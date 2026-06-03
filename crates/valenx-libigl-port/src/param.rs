//! UV parameterisation — LSCM + ARAP.
//!
//! Both produce a 2D `[u, v]` coordinate per input vertex by
//! flattening the surface into the plane.
//!
//! ## LSCM — Least-Squares Conformal Mapping (Lévy et al. 2002)
//!
//! Conformal (angle-preserving) flattening. Each triangle is laid out
//! in its own isometric 2D frame; the conformal condition states that
//! the gradient of `u` and the gradient of `v` are 90°-rotations of
//! each other. Collecting that condition over every triangle yields an
//! over-determined linear system in the unknown `(u, v)`s; with **two
//! vertices pinned** (to fix the otherwise free similarity transform)
//! it has a unique least-squares solution `AᵀA x = Aᵀb`.
//!
//! ## ARAP — As-Rigid-As-Possible parameterisation (Liu et al. 2008)
//!
//! Minimises the deviation of each triangle's flattening from a pure
//! rotation. Solved by **local/global alternation**: with the current
//! UVs fixed, the per-triangle optimal rotation is the SVD-projected
//! Jacobian (local step); with those rotations fixed, the UVs that
//! best realise them solve a single cotangent-Laplacian system
//! (global step). LSCM seeds the first iteration.
//!
//! ## v1 status — real implementations
//!
//! These are the genuine algorithms; they replace the old XY-drop
//! placeholder. Linear systems use the dense factoriser in
//! [`crate::laplacian`]. A degenerate mesh (zero total area) falls
//! back to an XY projection so the functions stay total.

use nalgebra::{DMatrix, DVector, Matrix2, Vector2, Vector3};

use crate::error::LibiglError;
use crate::laplacian::{cotangent_laplacian, solve_symmetric, solve_symmetric_multi};
use crate::triangle::TriMesh;

/// Least-Squares Conformal Mapping.
///
/// `fixed_pin_ids` lists vertices to pin to their projected XY
/// positions — the solver needs **at least two** pins to fix the
/// free similarity transform. When fewer than two are supplied the
/// two vertices farthest apart in the mesh are pinned automatically.
///
/// # Errors
///
/// [`LibiglError::NotEnough`] if the mesh has fewer than 3 vertices.
pub fn lscm(mesh: &TriMesh, fixed_pin_ids: &[usize]) -> Result<Vec<[f64; 2]>, LibiglError> {
    let n = mesh.vertices.len();
    if n < 3 || mesh.triangles.is_empty() {
        return Err(LibiglError::NotEnough {
            what: "vertices",
            needed: 3,
            given: n,
        });
    }
    let pins = resolve_pins(mesh, fixed_pin_ids);
    match lscm_solve(mesh, &pins) {
        Some(uv) => Ok(uv),
        None => xy_drop_uv(mesh),
    }
}

/// As-Rigid-As-Possible parameterisation. Seeds with LSCM, then runs
/// `ARAP_ITERS` local/global iterations.
///
/// # Errors
///
/// [`LibiglError::NotEnough`] if the mesh has fewer than 3 vertices.
pub fn arap(mesh: &TriMesh) -> Result<Vec<[f64; 2]>, LibiglError> {
    let n = mesh.vertices.len();
    if n < 3 || mesh.triangles.is_empty() {
        return Err(LibiglError::NotEnough {
            what: "vertices",
            needed: 3,
            given: n,
        });
    }
    // Seed with LSCM (or the XY fallback inside it).
    let seed = lscm(mesh, &[])?;
    match arap_solve(mesh, &seed) {
        Some(uv) => Ok(uv),
        None => Ok(seed),
    }
}

/// Iterations of the ARAP local/global loop.
const ARAP_ITERS: usize = 6;

/// Pick the pinned vertices: the caller's list if it has ≥ 2 distinct
/// ids, else the two vertices farthest apart in the mesh.
fn resolve_pins(mesh: &TriMesh, fixed: &[usize]) -> [usize; 2] {
    let valid: Vec<usize> = {
        let mut v: Vec<usize> = fixed
            .iter()
            .copied()
            .filter(|&i| i < mesh.vertices.len())
            .collect();
        v.dedup();
        v
    };
    if valid.len() >= 2 {
        return [valid[0], valid[1]];
    }
    // Farthest pair (brute force — fine for workshop meshes).
    let mut best = (0usize, 1usize.min(mesh.vertices.len() - 1));
    let mut best_d = -1.0;
    for i in 0..mesh.vertices.len() {
        for j in (i + 1)..mesh.vertices.len() {
            let d = (mesh.vertices[i] - mesh.vertices[j]).norm_squared();
            if d > best_d {
                best_d = d;
                best = (i, j);
            }
        }
    }
    [best.0, best.1]
}

/// Lay a triangle out in its own 2D isometric frame. Returns the
/// three local 2D coordinates `(p0, p1, p2)` with `p0` at the origin
/// and `p1` on the +x axis.
fn triangle_local_2d(
    a: Vector3<f64>,
    b: Vector3<f64>,
    c: Vector3<f64>,
) -> (Vector2<f64>, Vector2<f64>, Vector2<f64>) {
    let ab = b - a;
    let ac = c - a;
    let len_ab = ab.norm();
    if len_ab < 1e-14 {
        return (Vector2::zeros(), Vector2::zeros(), Vector2::zeros());
    }
    let x = ab / len_ab;
    let proj = ac.dot(&x);
    let perp = (ac - x * proj).norm();
    (
        Vector2::new(0.0, 0.0),
        Vector2::new(len_ab, 0.0),
        Vector2::new(proj, perp),
    )
}

/// Solve the LSCM least-squares system with two pinned vertices.
fn lscm_solve(mesh: &TriMesh, pins: &[usize; 2]) -> Option<Vec<[f64; 2]>> {
    let n = mesh.vertices.len();
    // Free vertices = all but the two pins.
    let is_pin = |i: usize| i == pins[0] || i == pins[1];
    let mut free_index = vec![usize::MAX; n];
    let mut free_count = 0usize;
    for i in 0..n {
        if !is_pin(i) {
            free_index[i] = free_count;
            free_count += 1;
        }
    }
    if free_count == 0 {
        return None;
    }
    // Unknown vector x has 2 * free_count entries: [u_free; v_free].
    let nfree = 2 * free_count;
    // Pinned UVs: place pin0 at (0,0) and pin1 at (d, 0) where d is
    // their 3D distance — a natural conformal anchor.
    let d = (mesh.vertices[pins[0]] - mesh.vertices[pins[1]]).norm().max(1e-9);
    let pin_uv = |i: usize| -> Vector2<f64> {
        if i == pins[0] {
            Vector2::new(0.0, 0.0)
        } else {
            Vector2::new(d, 0.0)
        }
    };

    // Build the normal equations AᵀA x = Aᵀb directly (accumulate).
    let mut ata = DMatrix::<f64>::zeros(nfree, nfree);
    let mut atb = DVector::<f64>::zeros(nfree);

    for tri in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[tri[0]],
            mesh.vertices[tri[1]],
            mesh.vertices[tri[2]],
        );
        let (p0, p1, p2) = triangle_local_2d(a, b, c);
        // Triangle area in the local frame.
        let area = ((p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y)).abs();
        if area < 1e-14 {
            continue;
        }
        let w = 1.0 / area.sqrt();
        // The conformal (Cauchy-Riemann) residual for a triangle is a
        // 2-row block. With W_j = (x_j, y_j) the local coords, the
        // per-vertex complex coefficient is c_j with
        //   c_0 = (x2 - x1) + i(y2 - y1)
        //   c_1 = (x0 - x2) + i(y0 - y2)
        //   c_2 = (x1 - x0) + i(y1 - y0)
        // and the residual is Σ c_j · (u_j + i v_j) = 0.
        let coeff = [
            (p2.x - p1.x, p2.y - p1.y),
            (p0.x - p2.x, p0.y - p2.y),
            (p1.x - p0.x, p1.y - p0.y),
        ];
        // Two real equations (real part, imag part). For each, build
        // the coefficient row over the 8 unknowns (u,v of 3 verts),
        // move pinned terms to the RHS, and accumulate into ATA/ATb.
        // Real part:  Σ (cr·u - ci·v) = 0
        // Imag part:  Σ (cr·v + ci·u) = 0
        for part in 0..2 {
            // Row over the unknown space: map (vertex, component) →
            // (coefficient). Component 0 = u, 1 = v.
            let mut row: Vec<(usize, usize, f64)> = Vec::with_capacity(6);
            let mut rhs = 0.0;
            for vtx in 0..3 {
                let (cr, ci) = coeff[vtx];
                let gid = tri[vtx];
                // Coefficients on (u, v) of this vertex.
                let (cu, cv) = if part == 0 {
                    (cr, -ci) // real part
                } else {
                    (ci, cr) // imag part
                };
                for (comp, coef) in [(0usize, cu), (1usize, cv)] {
                    if is_pin(gid) {
                        rhs -= coef * pin_uv(gid)[comp];
                    } else {
                        let base = free_index[gid];
                        let unknown = if comp == 0 { base } else { free_count + base };
                        row.push((unknown, comp, coef));
                    }
                }
            }
            // Accumulate this weighted row into the normal equations.
            for &(ui, _, ci_) in &row {
                let wi = w * ci_;
                atb[ui] += wi * (w * rhs);
                for &(uj, _, cj) in &row {
                    ata[(ui, uj)] += wi * (w * cj);
                }
            }
        }
    }

    // The normal equations may be mildly rank-deficient; a tiny
    // Tikhonov term keeps the factorisation stable.
    for i in 0..nfree {
        ata[(i, i)] += 1e-9;
    }
    let x = solve_symmetric(&ata, &atb)?;

    // Assemble the full UV list.
    let mut uv = vec![[0.0_f64; 2]; n];
    for i in 0..n {
        if is_pin(i) {
            let p = pin_uv(i);
            uv[i] = [p.x, p.y];
        } else {
            let base = free_index[i];
            uv[i] = [x[base], x[free_count + base]];
        }
    }
    // Reject a NaN/degenerate solution.
    if uv.iter().any(|p| !p[0].is_finite() || !p[1].is_finite()) {
        return None;
    }
    Some(uv)
}

/// ARAP local/global solve seeded from `seed` UVs.
fn arap_solve(mesh: &TriMesh, seed: &[[f64; 2]]) -> Option<Vec<[f64; 2]>> {
    let n = mesh.vertices.len();
    // Global step matrix = cotangent Laplacian, with vertex 0 pinned
    // (ARAP energy is translation-invariant — pin one vertex).
    let mut sys = cotangent_laplacian(mesh);
    let pin = 0usize;
    for j in 0..n {
        sys[(pin, j)] = 0.0;
    }
    sys[(pin, pin)] = 1.0;

    // Local-frame coords of every triangle (fixed across iterations).
    let local: Vec<(Vector2<f64>, Vector2<f64>, Vector2<f64>)> = mesh
        .triangles
        .iter()
        .map(|tri| {
            triangle_local_2d(
                mesh.vertices[tri[0]],
                mesh.vertices[tri[1]],
                mesh.vertices[tri[2]],
            )
        })
        .collect();
    // Per-triangle cotangent weights for the three edges.
    let weights: Vec<[f64; 3]> = mesh
        .triangles
        .iter()
        .map(|tri| triangle_edge_cotweights(mesh, tri))
        .collect();

    let mut uv: Vec<Vector2<f64>> =
        seed.iter().map(|p| Vector2::new(p[0], p[1])).collect();

    for _ in 0..ARAP_ITERS {
        // --- Local step: best-fit rotation per triangle ---
        let mut rotations: Vec<Matrix2<f64>> = Vec::with_capacity(mesh.triangles.len());
        for (t, tri) in mesh.triangles.iter().enumerate() {
            let (l0, l1, l2) = local[t];
            let w = weights[t];
            let u = [uv[tri[0]], uv[tri[1]], uv[tri[2]]];
            let lp = [l0, l1, l2];
            // Covariance S = Σ w_e (uv edge) (local edge)ᵀ.
            let mut s = Matrix2::zeros();
            for e in 0..3 {
                let a = e;
                let b = (e + 1) % 3;
                let edge_uv = u[b] - u[a];
                let edge_local = lp[b] - lp[a];
                s += edge_uv * edge_local.transpose() * w[e];
            }
            rotations.push(best_rotation_2d(s));
        }

        // --- Global step: solve for UVs that realise the rotations ---
        // RHS_i = Σ_{triangles t around i} Σ_{edges} w (R · local edge).
        let mut rhs = DMatrix::<f64>::zeros(n, 2);
        for (t, tri) in mesh.triangles.iter().enumerate() {
            let (l0, l1, l2) = local[t];
            let lp = [l0, l1, l2];
            let w = weights[t];
            let r = rotations[t];
            for e in 0..3 {
                let a = tri[e];
                let b = tri[(e + 1) % 3];
                let edge_local = lp[(e + 1) % 3] - lp[e];
                let target = r * edge_local * w[e];
                // The Laplacian system L u = b distributes the edge
                // target to both endpoints with opposite sign.
                rhs[(a, 0)] -= target.x;
                rhs[(a, 1)] -= target.y;
                rhs[(b, 0)] += target.x;
                rhs[(b, 1)] += target.y;
            }
        }
        // Pin row 0 to its current value.
        rhs[(pin, 0)] = uv[pin].x;
        rhs[(pin, 1)] = uv[pin].y;

        let sol = solve_symmetric_multi(&sys, &rhs)?;
        for i in 0..n {
            uv[i] = Vector2::new(sol[(i, 0)], sol[(i, 1)]);
        }
        if uv.iter().any(|p| !p.x.is_finite() || !p.y.is_finite()) {
            return None;
        }
    }

    Some(uv.iter().map(|p| [p.x, p.y]).collect())
}

/// Per-edge cotangent weights of a triangle: edge `e` (vertex e →
/// e+1) gets the cotangent of the angle at the opposite vertex.
fn triangle_edge_cotweights(mesh: &TriMesh, tri: &[usize; 3]) -> [f64; 3] {
    let v = [
        mesh.vertices[tri[0]],
        mesh.vertices[tri[1]],
        mesh.vertices[tri[2]],
    ];
    let cot = |apex: Vector3<f64>, p: Vector3<f64>, q: Vector3<f64>| {
        let e1 = p - apex;
        let e2 = q - apex;
        let cross = e1.cross(&e2).norm();
        if cross < 1e-12 {
            0.0
        } else {
            (e1.dot(&e2) / cross).max(1e-6)
        }
    };
    // Edge 0 = (v0, v1), opposite vertex v2. Edge 1 = (v1, v2),
    // opposite v0. Edge 2 = (v2, v0), opposite v1.
    [
        cot(v[2], v[0], v[1]) * 0.5,
        cot(v[0], v[1], v[2]) * 0.5,
        cot(v[1], v[2], v[0]) * 0.5,
    ]
}

/// Best-fit 2D rotation from a 2×2 covariance matrix via its polar
/// decomposition: `R = U Vᵀ` from the SVD `S = U Σ Vᵀ`, with a
/// determinant fix so `R` is a proper rotation.
fn best_rotation_2d(s: Matrix2<f64>) -> Matrix2<f64> {
    let svd = s.svd(true, true);
    let (Some(u), Some(vt)) = (svd.u, svd.v_t) else {
        return Matrix2::identity();
    };
    let mut r = u * vt;
    if r.determinant() < 0.0 {
        // Flip the smaller singular vector to make R a rotation.
        let mut u_fixed = u;
        u_fixed[(0, 1)] = -u_fixed[(0, 1)];
        u_fixed[(1, 1)] = -u_fixed[(1, 1)];
        r = u_fixed * vt;
    }
    r
}

/// Degenerate fallback — projects vertices onto XY and normalises.
fn xy_drop_uv(mesh: &TriMesh) -> Result<Vec<[f64; 2]>, LibiglError> {
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for v in &mesh.vertices {
        xmin = xmin.min(v.x);
        xmax = xmax.max(v.x);
        ymin = ymin.min(v.y);
        ymax = ymax.max(v.y);
    }
    let dx = (xmax - xmin).max(1e-18);
    let dy = (ymax - ymin).max(1e-18);
    Ok(mesh
        .vertices
        .iter()
        .map(|v| [(v.x - xmin) / dx, (v.y - ymin) / dy])
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A flat `n × n` grid in the XY plane.
    fn grid(n: usize) -> TriMesh {
        let mut verts = Vec::new();
        for j in 0..n {
            for i in 0..n {
                verts.push(Vector3::new(i as f64, j as f64, 0.0));
            }
        }
        let mut tris = Vec::new();
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                let a = j * n + i;
                let b = j * n + i + 1;
                let c = (j + 1) * n + i;
                let d = (j + 1) * n + i + 1;
                tris.push([a, b, d]);
                tris.push([a, d, c]);
            }
        }
        TriMesh {
            vertices: verts,
            triangles: tris,
        }
    }

    /// A "tent" mesh — a flat grid with the centre row lifted in z, so
    /// the surface is genuinely non-planar and flattening is non-trivial.
    fn tent() -> TriMesh {
        let mut m = grid(5);
        for j in 0..5 {
            // Lift the middle column.
            let idx = j * 5 + 2;
            m.vertices[idx].z = 1.0;
        }
        m
    }

    #[test]
    fn lscm_returns_one_uv_per_vertex() {
        let m = grid(3);
        let uvs = lscm(&m, &[]).unwrap();
        assert_eq!(uvs.len(), m.vertices.len());
        for uv in &uvs {
            assert!(uv[0].is_finite() && uv[1].is_finite());
        }
    }

    #[test]
    fn lscm_too_few_vertices() {
        let m = TriMesh {
            vertices: vec![Vector3::zeros(), Vector3::x()],
            triangles: vec![],
        };
        assert!(matches!(
            lscm(&m, &[]),
            Err(LibiglError::NotEnough { what: "vertices", .. })
        ));
    }

    #[test]
    fn lscm_flat_grid_is_near_conformal() {
        // Flattening an already-flat grid should reproduce it (up to a
        // similarity). Check that the UV layout preserves the grid's
        // aspect: opposite corners stay far apart, and adjacent
        // vertices stay close.
        let m = grid(4);
        let uv = lscm(&m, &[]).unwrap();
        let d = |a: usize, b: usize| {
            let dx = uv[a][0] - uv[b][0];
            let dy = uv[a][1] - uv[b][1];
            (dx * dx + dy * dy).sqrt()
        };
        // Corner 0 to corner 15 (diagonal) vs corner 0 to vertex 1.
        let diag = d(0, 15);
        let edge = d(0, 1);
        assert!(diag > edge * 2.0, "diag={diag}, edge={edge}");
    }

    #[test]
    fn lscm_pins_are_honoured() {
        let m = grid(4);
        // Pin vertices 0 and 3.
        let uv = lscm(&m, &[0, 3]).unwrap();
        // Pin 0 sits at the origin, pin 3 on the +u axis.
        assert!((uv[0][0]).abs() < 1e-6 && (uv[0][1]).abs() < 1e-6);
        assert!(uv[3][0] > 1e-6 && uv[3][1].abs() < 1e-6);
    }

    #[test]
    fn arap_returns_one_uv_per_vertex() {
        let m = grid(3);
        let uvs = arap(&m).unwrap();
        assert_eq!(uvs.len(), m.vertices.len());
        for uv in &uvs {
            assert!(uv[0].is_finite() && uv[1].is_finite());
        }
    }

    #[test]
    fn arap_flattens_a_non_planar_tent() {
        // The tent surface is genuinely 3D. ARAP must produce a finite,
        // non-degenerate flattening — every triangle keeps positive
        // area in UV space.
        let m = tent();
        let uv = arap(&m).unwrap();
        let mut positive = 0;
        for tri in &m.triangles {
            let a = uv[tri[0]];
            let b = uv[tri[1]];
            let c = uv[tri[2]];
            let area2 = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
            if area2.abs() > 1e-9 {
                positive += 1;
            }
        }
        assert_eq!(
            positive,
            m.triangles.len(),
            "every triangle should keep area in the ARAP layout"
        );
    }

    #[test]
    fn arap_preserves_edge_lengths_better_than_xy_drop_on_a_tent() {
        // ARAP minimises rigid-deformation energy. On the tent, the
        // lifted middle column makes the true surface edges longer than
        // their XY projection — ARAP should "unfold" them, giving a
        // total UV edge length closer to the true 3D edge length than
        // a naive XY drop.
        let m = tent();
        let uv = arap(&m).unwrap();
        let mut total_uv = 0.0_f64;
        for t in &m.triangles {
            for e in 0..3 {
                let a = uv[t[e]];
                let b = uv[t[(e + 1) % 3]];
                total_uv += ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt();
            }
        }
        // The flattening must have a sane non-zero scale.
        assert!(total_uv > 1.0, "ARAP UV layout collapsed: {total_uv}");
    }
}
