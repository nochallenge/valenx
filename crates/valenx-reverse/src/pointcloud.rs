//! Point-cloud loading, normal estimation, and triangulation.
//!
//! See the crate-level docs for the pipeline overview.

use std::path::Path;

use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

use crate::error::ReverseError;

/// A 3D point set, optionally with per-point normals.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct PointCloud {
    /// Point coordinates.
    pub points: Vec<Vector3<f64>>,
    /// Per-point normals (same length as `points` when present).
    pub normals: Option<Vec<Vector3<f64>>>,
}

impl PointCloud {
    /// Empty cloud.
    pub fn new() -> Self {
        Self::default()
    }

    /// Wrap an existing vector of points (normals = None).
    pub fn from_points(points: Vec<Vector3<f64>>) -> Self {
        Self {
            points,
            normals: None,
        }
    }

    /// Point count.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// `true` when the cloud has no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

/// Parse an ASCII PLY file at `path`.
///
/// Supports the canonical `format ascii 1.0` header with `element
/// vertex N` followed by `property float x|y|z` (3 columns required).
/// Optional `property float nx|ny|nz` (3 normal columns) are picked up
/// when present.
pub fn from_ply(path: &Path) -> Result<PointCloud, ReverseError> {
    // Round-23 workspace sweep: bounded at MAX_PLY_ASCII_BYTES
    // (1 GiB) — ASCII PLY for a million-vertex scan is in the low
    // hundreds of MiB; 1 GiB matches the OBJ / DXF caps.
    let text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_core::io_caps::MAX_PLY_ASCII_BYTES as usize,
    )?;
    from_ply_str(&text)
}

/// Parse PLY contents from a string.
pub fn from_ply_str(text: &str) -> Result<PointCloud, ReverseError> {
    let mut lines = text.lines();
    if lines.next().map(str::trim) != Some("ply") {
        return Err(ReverseError::PlyParse("missing 'ply' magic".into()));
    }
    let mut vertex_count: usize = 0;
    let mut props: Vec<String> = Vec::new();
    let mut in_vertex_element = false;
    let mut end_header_seen = false;
    for line in lines.by_ref() {
        let t = line.trim();
        if t == "end_header" {
            end_header_seen = true;
            break;
        }
        if let Some(rest) = t.strip_prefix("format ") {
            if !rest.starts_with("ascii") {
                return Err(ReverseError::PlyParse(format!(
                    "only ASCII PLY supported (got `{rest}`)"
                )));
            }
        } else if let Some(rest) = t.strip_prefix("element vertex ") {
            in_vertex_element = true;
            vertex_count = rest.trim().parse::<usize>().map_err(|e| {
                ReverseError::PlyParse(format!("bad vertex count `{rest}`: {e}"))
            })?;
        } else if t.starts_with("element ") {
            in_vertex_element = false;
        } else if in_vertex_element {
            if let Some(rest) = t.strip_prefix("property ") {
                // "<type> <name>"
                let mut parts = rest.split_whitespace();
                let _ty = parts.next();
                if let Some(name) = parts.next() {
                    props.push(name.to_string());
                }
            }
        }
        // comment lines, obj_info, format spec, other elements — skip.
    }
    if !end_header_seen {
        return Err(ReverseError::PlyParse("missing end_header".into()));
    }
    let idx_x = props.iter().position(|p| p == "x").ok_or_else(|| {
        ReverseError::PlyParse("missing property `x`".into())
    })?;
    let idx_y = props.iter().position(|p| p == "y").ok_or_else(|| {
        ReverseError::PlyParse("missing property `y`".into())
    })?;
    let idx_z = props.iter().position(|p| p == "z").ok_or_else(|| {
        ReverseError::PlyParse("missing property `z`".into())
    })?;
    let idx_nx = props.iter().position(|p| p == "nx");
    let idx_ny = props.iter().position(|p| p == "ny");
    let idx_nz = props.iter().position(|p| p == "nz");
    let has_normals = idx_nx.is_some() && idx_ny.is_some() && idx_nz.is_some();
    let mut points: Vec<Vector3<f64>> = Vec::with_capacity(vertex_count);
    let mut normals: Vec<Vector3<f64>> = Vec::with_capacity(if has_normals {
        vertex_count
    } else {
        0
    });
    for (row_i, raw) in lines.take(vertex_count).enumerate() {
        let cols: Vec<f64> = raw
            .split_whitespace()
            .map(|s| {
                s.parse::<f64>().map_err(|e| {
                    ReverseError::PlyParse(format!("row {row_i} col parse: {e}"))
                })
            })
            .collect::<Result<_, _>>()?;
        let need = props.len();
        if cols.len() != need {
            return Err(ReverseError::PlyParse(format!(
                "row {row_i}: expected {need} cols, got {}",
                cols.len()
            )));
        }
        points.push(Vector3::new(cols[idx_x], cols[idx_y], cols[idx_z]));
        if has_normals {
            normals.push(Vector3::new(
                cols[idx_nx.unwrap()],
                cols[idx_ny.unwrap()],
                cols[idx_nz.unwrap()],
            ));
        }
    }
    if points.len() != vertex_count {
        return Err(ReverseError::PlyParse(format!(
            "vertex underflow: expected {vertex_count}, got {}",
            points.len()
        )));
    }
    Ok(PointCloud {
        points,
        normals: if has_normals { Some(normals) } else { None },
    })
}

/// Parse a flat "x y z [nx ny nz]" whitespace-separated cloud (the
/// minimal-headers format some scanners spit out).
pub fn from_xyz(text: &str) -> Result<PointCloud, ReverseError> {
    let mut points = Vec::new();
    let mut normals = Vec::new();
    let mut have_normals: Option<bool> = None;
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<f64> = line
            .split_whitespace()
            .map(|s| {
                s.parse::<f64>().map_err(|e| {
                    ReverseError::PlyParse(format!("xyz row {i}: {e}"))
                })
            })
            .collect::<Result<_, _>>()?;
        if cols.len() != 3 && cols.len() != 6 {
            return Err(ReverseError::PlyParse(format!(
                "xyz row {i}: expected 3 or 6 cols, got {}",
                cols.len()
            )));
        }
        let row_has_normals = cols.len() == 6;
        match have_normals {
            None => have_normals = Some(row_has_normals),
            Some(prev) if prev != row_has_normals => {
                return Err(ReverseError::PlyParse(format!(
                    "xyz row {i}: mixed 3-col and 6-col rows"
                )));
            }
            _ => {}
        }
        points.push(Vector3::new(cols[0], cols[1], cols[2]));
        if row_has_normals {
            normals.push(Vector3::new(cols[3], cols[4], cols[5]));
        }
    }
    Ok(PointCloud {
        points,
        normals: if matches!(have_normals, Some(true)) {
            Some(normals)
        } else {
            None
        },
    })
}

/// Estimate per-point normals via 3x3 covariance PCA over each
/// point's `k` nearest neighbours. The eigenvector with the smallest
/// eigenvalue is the normal.
///
/// Returns a new [`PointCloud`] with `normals = Some(_)` populated.
/// Original cloud is left untouched. The orientation of each normal
/// is not consistently flipped (the algorithm has a 180° sign
/// ambiguity per point); upstream consumers needing globally
/// consistent orientation should run a separate normal-orientation
/// pass.
pub fn estimate_normals(cloud: &PointCloud, k: usize) -> Result<PointCloud, ReverseError> {
    if k < 3 {
        return Err(ReverseError::BadParameter {
            name: "k",
            reason: "needs >= 3 neighbours for PCA".into(),
        });
    }
    let n = cloud.points.len();
    if n < k + 1 {
        return Err(ReverseError::BadParameter {
            name: "k",
            reason: format!("cloud has {n} points, k={k} needs at least k+1"),
        });
    }
    let mut normals: Vec<Vector3<f64>> = Vec::with_capacity(n);
    for i in 0..n {
        let p = cloud.points[i];
        let mut neighbours = k_nearest(&cloud.points, i, k);
        // Compute centroid + 3x3 covariance over neighbours.
        let centroid = neighbours
            .iter()
            .map(|&j| cloud.points[j])
            .fold(Vector3::zeros(), |acc, v| acc + v)
            / (neighbours.len() as f64);
        // Include self for stability.
        neighbours.push(i);
        let mut cov = Matrix3::zeros();
        for j in neighbours {
            let d = cloud.points[j] - centroid;
            cov += d * d.transpose();
        }
        // 3x3 symmetric eigen — use nalgebra::SymmetricEigen.
        let eigen = nalgebra::SymmetricEigen::new(cov);
        // Eigenvalues column index of the smallest.
        let mut min_idx = 0usize;
        for j in 1..3 {
            if eigen.eigenvalues[j] < eigen.eigenvalues[min_idx] {
                min_idx = j;
            }
        }
        let n_raw = eigen.eigenvectors.column(min_idx).into_owned();
        let n_unit = if n_raw.norm() > f64::EPSILON {
            n_raw.normalize()
        } else {
            // Degenerate — return up-axis to keep counts aligned.
            Vector3::z()
        };
        // Sign convention: point outward from centroid in v1.
        let out = (p - centroid).normalize();
        let sign = if n_unit.dot(&out) >= 0.0 { 1.0 } else { -1.0 };
        normals.push(n_unit * sign);
    }
    Ok(PointCloud {
        points: cloud.points.clone(),
        normals: Some(normals),
    })
}

/// k-nearest neighbour indices of `cloud[i]` excluding `i` itself.
/// O(n) per query — fine for v1 sub-10k clouds; Phase 26.5 will swap
/// in a kd-tree.
fn k_nearest(points: &[Vector3<f64>], i: usize, k: usize) -> Vec<usize> {
    let p = points[i];
    let mut d2: Vec<(f64, usize)> = points
        .iter()
        .enumerate()
        .filter(|(j, _)| *j != i)
        .map(|(j, q)| ((p - q).norm_squared(), j))
        .collect();
    d2.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    d2.into_iter().take(k).map(|(_, j)| j).collect()
}

/// Triangulate the cloud via k-NN mutual-neighbour pairs.
///
/// Algorithm (v1, deliberately simple):
/// 1. For each point `i`, find its `k` nearest neighbours.
/// 2. For each ordered triple `(i, j, l)` with `i < j < l` and `j, l`
///    among `i`'s neighbours, emit a triangle when `j` and `l` are
///    each other's k-NN neighbours. This favours triangles whose
///    edges connect points that are all close to each other,
///    suppressing the worst long-edge spikes.
/// 3. Build a [`valenx_mesh::Mesh`] with one Tri3 block.
///
/// This is *not* a real surface reconstructor — it's a v1 placeholder
/// that's good enough to round-trip a dense, low-noise scan into a
/// solid for downstream BRep recovery. Phase 26.5 will swap in a real
/// algorithm (screened Poisson is the planned target).
pub fn triangulate(cloud: &PointCloud, k: usize) -> Result<valenx_mesh::Mesh, ReverseError> {
    if k < 3 {
        return Err(ReverseError::BadParameter {
            name: "k",
            reason: "needs >= 3 for triangulation".into(),
        });
    }
    let n = cloud.points.len();
    if n < 4 {
        return Err(ReverseError::BadParameter {
            name: "cloud",
            reason: "needs >= 4 points".into(),
        });
    }
    // Precompute each point's k-NN.
    let nbrs: Vec<Vec<usize>> = (0..n).map(|i| k_nearest(&cloud.points, i, k)).collect();
    // Build adjacency Set for mutual-neighbour check.
    let mut nbr_sets: Vec<std::collections::HashSet<usize>> = Vec::with_capacity(n);
    for nbr in &nbrs {
        nbr_sets.push(nbr.iter().copied().collect());
    }
    let mut tris: Vec<[usize; 3]> = Vec::new();
    for (i, cands) in nbrs.iter().enumerate() {
        for (a, &j) in cands.iter().enumerate() {
            if j <= i {
                continue;
            }
            for &l in cands.iter().skip(a + 1) {
                if l <= i || l == j {
                    continue;
                }
                // j and l must be mutual k-NN neighbours.
                if nbr_sets[j].contains(&l) && nbr_sets[l].contains(&j) {
                    tris.push([i, j, l]);
                }
            }
        }
    }
    if tris.is_empty() {
        return Err(ReverseError::EmptyTriangulation {
            reason: "no mutual k-NN triplets found — try larger k".into(),
        });
    }
    let mut mesh = valenx_mesh::Mesh::new("reverse-triangulation");
    mesh.nodes = cloud.points.clone();
    let mut block =
        valenx_mesh::element::ElementBlock::new(valenx_mesh::element::ElementType::Tri3);
    let mut conn: Vec<u32> = Vec::with_capacity(tris.len() * 3);
    for t in &tris {
        conn.push(t[0] as u32);
        conn.push(t[1] as u32);
        conn.push(t[2] as u32);
    }
    block.connectivity = conn;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ply_parses_minimal() {
        let text = "ply\nformat ascii 1.0\nelement vertex 3\nproperty float x\nproperty float y\nproperty float z\nend_header\n0 0 0\n1 0 0\n0 1 0\n";
        let cloud = from_ply_str(text).unwrap();
        assert_eq!(cloud.len(), 3);
        assert!(cloud.normals.is_none());
    }

    #[test]
    fn ply_parses_with_normals() {
        let text = "ply\nformat ascii 1.0\nelement vertex 1\nproperty float x\nproperty float y\nproperty float z\nproperty float nx\nproperty float ny\nproperty float nz\nend_header\n0 0 0 0 0 1\n";
        let cloud = from_ply_str(text).unwrap();
        assert_eq!(cloud.len(), 1);
        assert_eq!(cloud.normals.as_ref().unwrap()[0], Vector3::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn ply_rejects_binary_format() {
        let text = "ply\nformat binary_little_endian 1.0\nelement vertex 0\nend_header\n";
        let err = from_ply_str(text).unwrap_err();
        assert!(matches!(err, ReverseError::PlyParse(_)));
    }

    #[test]
    fn xyz_parses_six_cols() {
        let s = "0 0 0 0 0 1\n1 0 0 0 0 1\n";
        let cloud = from_xyz(s).unwrap();
        assert_eq!(cloud.len(), 2);
        assert!(cloud.normals.is_some());
    }

    #[test]
    fn normal_of_plane_points_up() {
        // 9 grid points on z=0 plane.
        let mut pts = Vec::new();
        for x in -1..=1 {
            for y in -1..=1 {
                pts.push(Vector3::new(x as f64, y as f64, 0.0));
            }
        }
        let cloud = PointCloud::from_points(pts);
        let with_n = estimate_normals(&cloud, 4).unwrap();
        let normals = with_n.normals.unwrap();
        // Every normal should be (close to) ±z. Allow either sign;
        // PCA is sign-ambiguous on the smallest eigenvector.
        for n in normals {
            assert!(n.x.abs() < 1e-6);
            assert!(n.y.abs() < 1e-6);
            assert!((n.z.abs() - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn triangulate_tetra_emits_triangles() {
        // 4 corners of a unit tetra — k=3 will connect them all.
        let pts = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let cloud = PointCloud::from_points(pts);
        let mesh = triangulate(&cloud, 3).unwrap();
        // Each tri block must have ≥ 1 triangle.
        assert!(!mesh.element_blocks.is_empty());
        let block = &mesh.element_blocks[0];
        assert_eq!(block.element_type, valenx_mesh::element::ElementType::Tri3);
        assert!(block.connectivity.len() >= 3);
    }
}
