//! Pre-flight planner: validate inputs + build grids, but skip search.

use serde::Serialize;

use crate::config::DockConfig;
use crate::error::DockError;
use crate::grid::GridBundle;
use crate::ligand::Ligand;
use crate::receptor::Receptor;

/// Result of [`dock_dry_run`]: enough to let a caller decide whether
/// the configured job is sensible before actually running it.
#[derive(Debug, Clone, Serialize)]
pub struct DryRunPlan {
    /// Receptor atom count.
    pub receptor_atoms: usize,
    /// Ligand atom count (all, including H).
    pub ligand_atoms: usize,
    /// Ligand heavy-atom count.
    pub ligand_heavy_atoms: usize,
    /// Number of rotatable bonds.
    pub n_torsions: usize,
    /// Grid dimensions (nx, ny, nz).
    pub grid_dims: (usize, usize, usize),
    /// Total grid voxel count across all atom-type maps.
    pub total_voxels: usize,
    /// Rough estimate of search wall-clock seconds, derived from
    /// exhaustiveness × n_torsions on a calibration baseline.
    pub estimated_seconds: f64,
    /// Notes the caller may want to surface (out-of-box ligand, etc.).
    pub warnings: Vec<String>,
}

/// Validate inputs, precompute grids, and report an execution plan
/// without running the search.
pub fn dock_dry_run(
    receptor_pdbqt: &str,
    ligand_pdbqt: &str,
    config: &DockConfig,
) -> Result<DryRunPlan, DockError> {
    use crate::atom_type::Ad4AtomType;
    config.validate()?;
    let receptor = Receptor::from_pdbqt(receptor_pdbqt)?;
    let ligand = Ligand::from_pdbqt(ligand_pdbqt)?;
    let heavy = ligand
        .atoms
        .iter()
        .filter(|a| !matches!(a.ad4_type, Ad4AtomType::H | Ad4AtomType::HD))
        .count();
    let grids = GridBundle::build(
        &receptor,
        &ligand,
        config.grid_origin(),
        config.grid_spacing,
        config.grid_dims(),
    );
    let total_voxels: usize = grids
        .grids
        .values()
        .map(|g| g.dims.0 * g.dims.1 * g.dims.2)
        .sum();
    // Empirical baseline: ~1 sec per (exhaustiveness × max(1, n_torsions))
    // on a modern laptop with release builds.
    let estimated_seconds =
        config.exhaustiveness as f64 * (ligand.n_torsions().max(1)) as f64 * 1.0;
    let mut warnings = Vec::new();
    if heavy == 0 {
        warnings.push("ligand has no heavy atoms".into());
    }
    if receptor.atoms.is_empty() {
        warnings.push("receptor parsed empty".into());
    }
    let lig_centroid = ligand
        .atoms
        .iter()
        .fold(nalgebra::Vector3::zeros(), |s, a| s + a.position)
        / ligand.atoms.len().max(1) as f64;
    let half = config.size / 2.0;
    if (lig_centroid.x - config.center.x).abs() > half.x
        || (lig_centroid.y - config.center.y).abs() > half.y
        || (lig_centroid.z - config.center.z).abs() > half.z
    {
        warnings.push("ligand initial centroid is outside the search box (search will still run; the random restarts seed inside the box)".into());
    }
    Ok(DryRunPlan {
        receptor_atoms: receptor.atoms.len(),
        ligand_atoms: ligand.atoms.len(),
        ligand_heavy_atoms: heavy,
        n_torsions: ligand.n_torsions(),
        grid_dims: config.grid_dims(),
        total_voxels,
        estimated_seconds,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    // PDBQT requires lines >= 79 chars (ad4_type lives in cols 78..79).
    // Build fixtures programmatically so editor/markdown trimming can't
    // strip the trailing whitespace and break the parser.
    fn atom_line(serial: i32, name: &str, x: f64, y: f64, z: f64, ad4: &str) -> String {
        // Field columns per PDBQT spec: 1-6 record, 7-11 serial,
        // 13-16 name, 18-20 residue, 22 chain, 23-26 resseq,
        // 31-38 x, 39-46 y, 47-54 z, 55-60 occ, 61-66 b, 71-76 q, 78-79 ad4.
        let line = format!(
            "ATOM  {serial:>5} {name:<4} LIG A{:>4}    {x:8.3}{y:8.3}{z:8.3}  1.00  0.00     0.000 {ad4:>2}",
            1,
        );
        // Sanity: must reach col 79.
        assert!(
            line.len() >= 79,
            "atom_line produced {} chars: {:?}",
            line.len(),
            line
        );
        line
    }

    #[test]
    fn dry_run_reports_atom_counts_and_voxels() {
        let receptor = format!(
            "{}\n{}\n",
            atom_line(1, "CA", 0.0, 0.0, 0.0, "C"),
            atom_line(2, "CB", 1.5, 0.0, 0.0, "C"),
        );
        let ligand = format!(
            "ROOT\n{}\n{}\nENDROOT\nTORSDOF 0\n",
            atom_line(1, "C1", 0.0, 0.0, 0.0, "C"),
            atom_line(2, "O1", 1.5, 0.0, 0.0, "OA"),
        );
        let cfg = DockConfig {
            center: Vector3::new(0.75, 0.0, 0.0),
            size: Vector3::new(5.0, 5.0, 5.0),
            ..DockConfig::default()
        };
        let plan = dock_dry_run(&receptor, &ligand, &cfg).unwrap();
        assert_eq!(plan.receptor_atoms, 2);
        assert_eq!(plan.ligand_atoms, 2);
        assert_eq!(plan.ligand_heavy_atoms, 2);
        assert_eq!(plan.n_torsions, 0);
        assert!(plan.total_voxels > 0);
    }

    #[test]
    fn dry_run_warns_when_ligand_outside_box() {
        let receptor = format!("{}\n", atom_line(1, "CA", 0.0, 0.0, 0.0, "C"));
        let ligand = format!(
            "ROOT\n{}\nENDROOT\nTORSDOF 0\n",
            atom_line(1, "C1", 100.0, 0.0, 0.0, "C"),
        );
        let cfg = DockConfig {
            center: Vector3::zeros(),
            size: Vector3::new(5.0, 5.0, 5.0),
            ..DockConfig::default()
        };
        let plan = dock_dry_run(&receptor, &ligand, &cfg).unwrap();
        assert!(plan
            .warnings
            .iter()
            .any(|w| w.contains("outside the search box")));
    }
}
