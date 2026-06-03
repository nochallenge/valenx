//! FEM result parsing + visualization tessellator.

use std::path::Path;

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::solver::FemSolverChoice;

/// Errors from the FEM post-processor.
#[derive(Debug, Error)]
pub enum FemResultError {
    /// IO failure reading a solver output file.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    /// Parse failure (malformed solver output).
    #[error("parse: {0}")]
    Parse(String),
    /// Solver-specific parsing path not yet wired.
    #[error("parser for {0:?} is not yet implemented")]
    UnsupportedSolver(FemSolverChoice),
}

/// Parsed FEM analysis result.
#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct FemResult {
    /// Per-node displacement field [ux, uy, uz] in metres.
    pub displacement_field: Vec<[f64; 3]>,
    /// Per-node scalar von Mises stress in Pa.
    pub stress_field: Vec<f64>,
    /// Tessellated mesh (typically identical to the input surface
    /// mesh — included so the post-processor is self-contained).
    pub mesh: Mesh,
}

impl FemResult {
    /// Parse solver output at `path` for the given solver. v1 reads
    /// CalculiX `.frd` (text) by grepping for displacement + stress
    /// blocks; Elmer / Code_Aster parsing returns
    /// [`FemResultError::UnsupportedSolver`] (the existing adapters
    /// own their own parsers — this stub leaves room for direct
    /// in-process parsing later).
    ///
    /// # Errors
    ///
    /// - [`FemResultError::Io`] for read failures.
    /// - [`FemResultError::Parse`] for malformed output.
    /// - [`FemResultError::UnsupportedSolver`] when the solver's
    ///   parser is not wired in v1.
    pub fn from_solver_output(
        path: &Path,
        solver_kind: FemSolverChoice,
    ) -> Result<Self, FemResultError> {
        match solver_kind {
            FemSolverChoice::CalculiX => parse_frd(path),
            other => Err(FemResultError::UnsupportedSolver(other)),
        }
    }

    /// Maximum displacement magnitude across the field. Useful for
    /// the scale slider's "auto" preset.
    pub fn max_displacement(&self) -> f64 {
        self.displacement_field
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0f64, f64::max)
    }

    /// Maximum stress value across the field.
    pub fn max_stress(&self) -> f64 {
        self.stress_field.iter().copied().fold(0.0f64, f64::max)
    }

    /// Tessellate this result as a coloured mesh — duplicate the
    /// stored mesh but colour each vertex by stress (returned as a
    /// `(mesh, per_vertex_normalised_stress)` pair so a renderer can
    /// look up a colour ramp).
    pub fn tessellate_visualization(&self) -> (Mesh, Vec<f64>) {
        let max = self.max_stress().max(1e-9);
        let normalised: Vec<f64> = self.stress_field.iter().map(|s| s / max).collect();
        (self.mesh.clone(), normalised)
    }

    /// Animation step — interpolate the mesh from undeformed (`t=0`)
    /// to deformed (`t=1`). Returns a fresh mesh whose nodes are
    /// `mesh.nodes[i] + t * displacement_field[i]`.
    pub fn animate(&self, t: f64) -> Mesh {
        let t = t.clamp(0.0, 1.0);
        let mut out = self.mesh.clone();
        for (i, node) in out.nodes.iter_mut().enumerate() {
            if let Some(d) = self.displacement_field.get(i) {
                node.x += t * d[0];
                node.y += t * d[1];
                node.z += t * d[2];
            }
        }
        out
    }
}

/// Parse a CalculiX `.frd` text file. v1 is a forgiving line scan:
/// look for the displacement (`-4 DISP`) block and the stress block
/// (`-4 STRESS`). The full spec is much larger; we extract the two
/// fields the FEM panel needs.
///
/// Round-18 M1: read is capped at
/// [`valenx_core::io_caps::MAX_FRD_FILE_BYTES`] to refuse a hostile
/// multi-GB `.frd` before allocation.
fn parse_frd(path: &Path) -> Result<FemResult, FemResultError> {
    let text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_core::io_caps::MAX_FRD_FILE_BYTES,
    )?;
    let mut displacement: Vec<[f64; 3]> = Vec::new();
    let mut stresses: Vec<f64> = Vec::new();
    let nodes: Vec<Vector3<f64>> = Vec::new();
    let connectivity: Vec<u32> = Vec::new();
    let mut in_disp = false;
    let mut in_stress = false;
    for line in text.lines() {
        if line.contains("-4 DISP") {
            in_disp = true;
            in_stress = false;
            continue;
        }
        if line.contains("-4 STRESS") {
            in_stress = true;
            in_disp = false;
            continue;
        }
        if line.starts_with(" -3") {
            in_disp = false;
            in_stress = false;
            continue;
        }
        if in_disp && line.starts_with(" -1") {
            // Format: " -1{node_id} {ux} {uy} {uz}"
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 5 {
                let ux = parts[2].parse().unwrap_or(0.0);
                let uy = parts[3].parse().unwrap_or(0.0);
                let uz = parts[4].parse().unwrap_or(0.0);
                displacement.push([ux, uy, uz]);
            }
        }
        if in_stress && line.starts_with(" -1") {
            // Take the first 6 stress components; combine into a
            // crude scalar von-Mises proxy.
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 8 {
                let sx: f64 = parts[2].parse().unwrap_or(0.0);
                let sy: f64 = parts[3].parse().unwrap_or(0.0);
                let sz: f64 = parts[4].parse().unwrap_or(0.0);
                let sxy: f64 = parts[5].parse().unwrap_or(0.0);
                let syz: f64 = parts[6].parse().unwrap_or(0.0);
                let szx: f64 = parts[7].parse().unwrap_or(0.0);
                // von Mises stress: sqrt(0.5 * ((sx-sy)^2 + (sy-sz)^2
                // + (sz-sx)^2) + 3*(sxy^2 + syz^2 + szx^2))
                let vm = (0.5 * ((sx - sy).powi(2) + (sy - sz).powi(2) + (sz - sx).powi(2))
                    + 3.0 * (sxy * sxy + syz * syz + szx * szx))
                    .sqrt();
                stresses.push(vm);
            }
        }
    }
    // For v1 we don't reconstruct the mesh — leave nodes / blocks
    // empty and let the caller pair this Result with the source
    // tessellation. The presence of `nodes` + `connectivity` here is
    // forward-compatible scaffolding.
    let mesh = if !nodes.is_empty() && !connectivity.is_empty() {
        let mut m = Mesh::new("frd_mesh");
        m.nodes = nodes;
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = connectivity;
        m.element_blocks.push(block);
        m
    } else {
        Mesh::new("frd_empty")
    };
    Ok(FemResult {
        displacement_field: displacement,
        stress_field: stresses,
        mesh,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_displacement_finds_largest_magnitude() {
        let r = FemResult {
            displacement_field: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 4.0, 0.0]],
            stress_field: vec![],
            mesh: Mesh::new("t"),
        };
        assert!((r.max_displacement() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn max_stress_finds_max() {
        let r = FemResult {
            displacement_field: vec![],
            stress_field: vec![1.0, 5.0, 2.0],
            mesh: Mesh::new("t"),
        };
        assert!((r.max_stress() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn tessellate_visualization_normalises_stress() {
        let r = FemResult {
            displacement_field: vec![],
            stress_field: vec![10.0, 5.0, 0.0],
            mesh: Mesh::new("t"),
        };
        let (_, normalised) = r.tessellate_visualization();
        assert!((normalised[0] - 1.0).abs() < 1e-9);
        assert!((normalised[1] - 0.5).abs() < 1e-9);
        assert!((normalised[2] - 0.0).abs() < 1e-9);
    }

    #[test]
    fn animate_at_zero_returns_undeformed() {
        let mut mesh = Mesh::new("t");
        mesh.nodes.push(Vector3::new(1.0, 2.0, 3.0));
        let r = FemResult {
            displacement_field: vec![[10.0, 20.0, 30.0]],
            stress_field: vec![],
            mesh,
        };
        let m0 = r.animate(0.0);
        assert!((m0.nodes[0].x - 1.0).abs() < 1e-9);
    }

    #[test]
    fn animate_at_one_returns_full_deformation() {
        let mut mesh = Mesh::new("t");
        mesh.nodes.push(Vector3::new(1.0, 2.0, 3.0));
        let r = FemResult {
            displacement_field: vec![[10.0, 20.0, 30.0]],
            stress_field: vec![],
            mesh,
        };
        let m1 = r.animate(1.0);
        assert!((m1.nodes[0].x - 11.0).abs() < 1e-9);
        assert!((m1.nodes[0].y - 22.0).abs() < 1e-9);
        assert!((m1.nodes[0].z - 33.0).abs() < 1e-9);
    }

    #[test]
    fn animate_clamps_t_to_unit_range() {
        let mut mesh = Mesh::new("t");
        mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        let r = FemResult {
            displacement_field: vec![[1.0, 0.0, 0.0]],
            stress_field: vec![],
            mesh,
        };
        let m_neg = r.animate(-1.0);
        let m_big = r.animate(2.0);
        assert!((m_neg.nodes[0].x - 0.0).abs() < 1e-9);
        assert!((m_big.nodes[0].x - 1.0).abs() < 1e-9);
    }

    #[test]
    fn unsupported_solver_returns_err() {
        let path = std::env::temp_dir().join("valenx_fem_test_does_not_exist");
        let _ = std::fs::remove_file(&path);
        // Code_Aster path not wired in v1.
        let err = FemResult::from_solver_output(&path, FemSolverChoice::CodeAster).unwrap_err();
        assert!(matches!(err, FemResultError::UnsupportedSolver(_)));
    }

    #[test]
    fn parse_frd_extracts_displacement() {
        // A raw string literal — the CalculiX FRD format is
        // column-sensitive: data records begin with a leading space then
        // the record key (" -1", " -3"), and `parse_frd` keys off
        // exactly that. A `\`-continued string literal would strip the
        // leading whitespace and the records would never be recognised.
        let frd = " -4 DISP\n -1   1 0.0 0.0 0.0\n -1   2 1.0 0.0 0.0\n -3\n";
        let tmp = std::env::temp_dir().join("valenx_fem_frd_test.frd");
        std::fs::write(&tmp, frd).unwrap();
        let r = parse_frd(&tmp).unwrap();
        let _ = std::fs::remove_file(&tmp);
        assert_eq!(r.displacement_field.len(), 2);
        assert!((r.displacement_field[1][0] - 1.0).abs() < 1e-9);
    }

    /// Round-18 M1 RED→GREEN: `parse_frd` must refuse a `.frd` larger
    /// than `MAX_FRD_FILE_BYTES`. Sparse-file trick — seek to past the
    /// cap and write one byte so `metadata.len()` reports the size
    /// without us actually writing 256 MiB+ of zeros.
    #[test]
    fn parse_frd_rejects_oversize_file() {
        use std::io::{Seek, SeekFrom, Write};
        let tmp = std::env::temp_dir().join(format!(
            "valenx_fem_frd_oversize-{}.frd",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&tmp)
                .unwrap();
            // 1 byte past the FRD cap → metadata.len() reports oversize.
            f.seek(SeekFrom::Start(
                valenx_core::io_caps::MAX_FRD_FILE_BYTES as u64 + 1,
            ))
            .unwrap();
            f.write_all(b"x").unwrap();
        }
        let err = parse_frd(&tmp).expect_err("must reject oversize .frd");
        let _ = std::fs::remove_file(&tmp);
        match err {
            FemResultError::Io(io_err) => {
                assert_eq!(
                    io_err.kind(),
                    std::io::ErrorKind::InvalidData,
                    "expected InvalidData (cap rejection), got: {io_err}"
                );
            }
            other => panic!("expected Io(InvalidData), got: {other:?}"),
        }
    }
}
