//! Feature 19 — binding-pocket detection.
//!
//! Before you can dock you need a *where* — the search box has to sit
//! on an actual cavity. fpocket finds cavities by packing "alpha
//! spheres" (spheres touching four protein atoms with no atom inside)
//! into the protein surface and clustering them.
//!
//! This module implements a geometric pocket detector in that spirit,
//! using a *grid probe* rather than the alpha-sphere Voronoi
//! construction:
//!
//! 1. lay a regular grid over the protein's bounding box;
//! 2. keep a grid point if it is **buried** — not inside any atom, but
//!    surrounded by protein atoms in enough directions to be a cavity
//!    rather than open solvent;
//! 3. cluster the kept points by spatial adjacency;
//! 4. each cluster is a [`Pocket`] — its centre, volume estimate, a
//!    druggability-style score and the lining residues.
//!
//! This finds real cavities and ranks them sensibly; it is *not* the
//! exact fpocket alpha-sphere algorithm (see the crate-level v1 note).

use std::collections::VecDeque;

use nalgebra::Vector3;
use valenx_biostruct::structure::Structure;

use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;

/// A detected binding pocket.
#[derive(Clone, Debug, PartialEq)]
pub struct Pocket {
    /// Geometric centre of the pocket (Å).
    pub center: Vector3<f64>,
    /// Estimated pocket volume (Å³) — the kept grid-point count times
    /// the per-point cell volume.
    pub volume: f64,
    /// Number of buried grid points that formed the pocket.
    pub point_count: usize,
    /// A heuristic druggability / quality score in `[0, 1]` — larger
    /// for deeper, more enclosed, better-sized pockets.
    pub score: f64,
    /// `(chain_id, residue_name, seq_num)` of every residue lining the
    /// pocket.
    pub lining_residues: Vec<(String, String, i32)>,
}

impl Pocket {
    /// A [`GridBox`] sized to enclose this pocket, with `padding` Å of
    /// margin — ready to hand to a docking run.
    pub fn to_grid_box(&self, padding: f64) -> Result<GridBox> {
        // Estimate an edge length from the volume (a cube of equal
        // volume), padded.
        let edge = self.volume.max(1.0).cbrt() + 2.0 * padding.max(0.0);
        GridBox::cubic([self.center.x, self.center.y, self.center.z], edge.max(6.0))
    }
}

/// Tunable parameters for the pocket detector.
#[derive(Clone, Copy, Debug)]
pub struct PocketParams {
    /// Grid spacing (Å). Finer grids find more detail but cost more.
    pub spacing: f64,
    /// A grid point counts as "buried" when at least this many of the
    /// 14 probe directions hit a protein atom within `probe_reach`.
    pub min_buried_directions: usize,
    /// How far (Å) a burial probe ray reaches looking for protein.
    pub probe_reach: f64,
    /// Discard clusters smaller than this many grid points (filters
    /// shallow surface dimples).
    pub min_cluster_points: usize,
}

impl Default for PocketParams {
    fn default() -> Self {
        PocketParams {
            spacing: 1.0,
            min_buried_directions: 9,
            probe_reach: 8.0,
            min_cluster_points: 12,
        }
    }
}

/// The 14 burial-probe directions: the 6 face and 8 corner directions
/// of a cube. A grid point surrounded in most of these is enclosed.
fn probe_directions() -> [Vector3<f64>; 14] {
    let s = 1.0 / 3.0_f64.sqrt();
    [
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, -1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 0.0, -1.0),
        Vector3::new(s, s, s),
        Vector3::new(s, s, -s),
        Vector3::new(s, -s, s),
        Vector3::new(s, -s, -s),
        Vector3::new(-s, s, s),
        Vector3::new(-s, s, -s),
        Vector3::new(-s, -s, s),
        Vector3::new(-s, -s, -s),
    ]
}

/// A heavy protein atom for pocket detection.
struct ProtAtom {
    pos: Vector3<f64>,
    radius: f64,
    chain: String,
    residue: String,
    seq_num: i32,
}

/// Approximate van der Waals radius from an element symbol.
fn element_radius(element: &str) -> f64 {
    match element.to_ascii_uppercase().as_str() {
        "C" => 1.7,
        "N" => 1.55,
        "O" => 1.52,
        "S" => 1.8,
        "P" => 1.8,
        "H" | "D" => 1.2,
        "FE" | "ZN" | "MG" | "CA" | "MN" | "CU" => 1.4,
        _ => 1.7,
    }
}

/// Feature 19 — detect binding pockets in a receptor structure.
///
/// Returns pockets sorted by score (best first). An empty result
/// means no buried cavity passed the `min_cluster_points` filter.
///
/// Returns [`DockScreenError::InvalidReceptor`] if the structure has
/// no heavy atoms.
pub fn detect_pockets(receptor: &Structure, params: &PocketParams) -> Result<Vec<Pocket>> {
    if !params.spacing.is_finite() || params.spacing <= 0.0 {
        return Err(DockScreenError::invalid("spacing", "must be positive"));
    }
    // Collect heavy protein atoms.
    let model = receptor.first_model();
    let mut atoms: Vec<ProtAtom> = Vec::new();
    for chain in &model.chains {
        for residue in &chain.residues {
            for atom in &residue.atoms {
                if atom.is_hydrogen() {
                    continue;
                }
                atoms.push(ProtAtom {
                    pos: Vector3::new(atom.coord.x, atom.coord.y, atom.coord.z),
                    radius: element_radius(&atom.element),
                    chain: chain.id.clone(),
                    residue: residue.name.clone(),
                    seq_num: residue.seq_num,
                });
            }
        }
    }
    if atoms.is_empty() {
        return Err(DockScreenError::invalid_receptor(
            "structure has no heavy atoms",
        ));
    }

    // Bounding box of the protein, padded by the probe reach.
    let mut lo = atoms[0].pos;
    let mut hi = atoms[0].pos;
    for a in &atoms {
        lo = lo.inf(&a.pos);
        hi = hi.sup(&a.pos);
    }
    let pad = Vector3::repeat(2.0);
    lo -= pad;
    hi += pad;

    let dirs = probe_directions();
    let sp = params.spacing;
    let nx = (((hi.x - lo.x) / sp).ceil() as usize).max(1);
    let ny = (((hi.y - lo.y) / sp).ceil() as usize).max(1);
    let nz = (((hi.z - lo.z) / sp).ceil() as usize).max(1);

    // For each grid point, decide whether it is a buried cavity point.
    let buried = |p: Vector3<f64>| -> bool {
        // Inside an atom → not a cavity point.
        for a in &atoms {
            if (p - a.pos).norm() < a.radius {
                return false;
            }
        }
        // Count probe directions that hit protein within reach.
        let mut hits = 0usize;
        for d in &dirs {
            let mut t = sp;
            let mut hit = false;
            while t <= params.probe_reach {
                let q = p + d * t;
                for a in &atoms {
                    if (q - a.pos).norm() < a.radius + 0.5 {
                        hit = true;
                        break;
                    }
                }
                if hit {
                    break;
                }
                t += sp;
            }
            if hit {
                hits += 1;
            }
        }
        hits >= params.min_buried_directions
    };

    // 3-D array of buried flags.
    let idx = |ix: usize, iy: usize, iz: usize| ix + iy * nx + iz * nx * ny;
    let mut flags = vec![false; nx * ny * nz];
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let p = lo + Vector3::new(ix as f64, iy as f64, iz as f64) * sp;
                if buried(p) {
                    flags[idx(ix, iy, iz)] = true;
                }
            }
        }
    }

    // Cluster buried points by 6-connectivity (flood fill).
    let mut visited = vec![false; nx * ny * nz];
    let cell_volume = sp * sp * sp;
    let mut pockets: Vec<Pocket> = Vec::new();
    for sz in 0..nz {
        for sy in 0..ny {
            for sx in 0..nx {
                let start = idx(sx, sy, sz);
                if !flags[start] || visited[start] {
                    continue;
                }
                // BFS over the connected buried region.
                let mut cluster: Vec<(usize, usize, usize)> = Vec::new();
                let mut queue: VecDeque<(usize, usize, usize)> = VecDeque::new();
                visited[start] = true;
                queue.push_back((sx, sy, sz));
                while let Some((cx, cy, cz)) = queue.pop_front() {
                    cluster.push((cx, cy, cz));
                    let mut nbrs: Vec<(isize, isize, isize)> = Vec::new();
                    nbrs.push((cx as isize + 1, cy as isize, cz as isize));
                    nbrs.push((cx as isize - 1, cy as isize, cz as isize));
                    nbrs.push((cx as isize, cy as isize + 1, cz as isize));
                    nbrs.push((cx as isize, cy as isize - 1, cz as isize));
                    nbrs.push((cx as isize, cy as isize, cz as isize + 1));
                    nbrs.push((cx as isize, cy as isize, cz as isize - 1));
                    for (mx, my, mz) in nbrs {
                        if mx < 0
                            || my < 0
                            || mz < 0
                            || mx as usize >= nx
                            || my as usize >= ny
                            || mz as usize >= nz
                        {
                            continue;
                        }
                        let ni = idx(mx as usize, my as usize, mz as usize);
                        if flags[ni] && !visited[ni] {
                            visited[ni] = true;
                            queue.push_back((mx as usize, my as usize, mz as usize));
                        }
                    }
                }
                if cluster.len() < params.min_cluster_points {
                    continue;
                }
                pockets.push(build_pocket(&cluster, lo, sp, cell_volume, &atoms));
            }
        }
    }
    pockets.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(pockets)
}

/// Assemble a [`Pocket`] from a cluster of buried grid cells.
fn build_pocket(
    cluster: &[(usize, usize, usize)],
    origin: Vector3<f64>,
    spacing: f64,
    cell_volume: f64,
    atoms: &[ProtAtom],
) -> Pocket {
    // Centre = mean of the cluster's cell positions.
    let mut centre = Vector3::zeros();
    for &(ix, iy, iz) in cluster {
        centre += origin + Vector3::new(ix as f64, iy as f64, iz as f64) * spacing;
    }
    centre /= cluster.len() as f64;
    let volume = cluster.len() as f64 * cell_volume;

    // Lining residues: every residue with a heavy atom within 5 Å of
    // any cluster cell.
    let mut lining: Vec<(String, String, i32)> = Vec::new();
    let lining_cut2 = 5.0 * 5.0;
    for a in atoms {
        let near = cluster.iter().any(|&(ix, iy, iz)| {
            let p = origin + Vector3::new(ix as f64, iy as f64, iz as f64) * spacing;
            (p - a.pos).norm_squared() <= lining_cut2
        });
        if near {
            let key = (a.chain.clone(), a.residue.clone(), a.seq_num);
            if !lining.contains(&key) {
                lining.push(key);
            }
        }
    }

    // Heuristic score: reward a sensible volume (a real druggable
    // pocket is ~200-1000 Å³) and a good lining-residue count.
    let volume_term = {
        // Triangular preference peaking around 600 Å³.
        let v = volume;
        if v <= 0.0 {
            0.0
        } else if v < 600.0 {
            v / 600.0
        } else {
            (1.0 - (v - 600.0) / 2400.0).max(0.0)
        }
    };
    let residue_term = (lining.len() as f64 / 25.0).min(1.0);
    let score = (0.6 * volume_term + 0.4 * residue_term).clamp(0.0, 1.0);

    Pocket {
        center: centre,
        volume,
        point_count: cluster.len(),
        score,
        lining_residues: lining,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_biostruct::load_structure;

    /// Build a PDB string for a hollow box of carbon atoms — the
    /// interior is a clean cavity the detector should find. The shell
    /// is a 7×7×7 lattice (343 atoms) with the centre 3×3×3 removed.
    fn hollow_box_pdb() -> String {
        let mut s = String::new();
        let mut serial = 1;
        for ix in 0..7 {
            for iy in 0..7 {
                for iz in 0..7 {
                    // Hollow out the central 3×3×3.
                    let inner =
                        (2..=4).contains(&ix) && (2..=4).contains(&iy) && (2..=4).contains(&iz);
                    if inner {
                        continue;
                    }
                    let x = ix as f64 * 3.0;
                    let y = iy as f64 * 3.0;
                    let z = iz as f64 * 3.0;
                    s.push_str(&format!(
                        "ATOM  {serial:>5}  CA  ALA A{serial:>4}    {x:8.3}{y:8.3}{z:8.3}  1.00  0.00           C\n"
                    ));
                    serial += 1;
                }
            }
        }
        s.push_str("END\n");
        s
    }

    #[test]
    fn rejects_bad_spacing() {
        let s = load_structure(&hollow_box_pdb(), "x").unwrap();
        let params = PocketParams {
            spacing: 0.0,
            ..PocketParams::default()
        };
        assert!(detect_pockets(&s, &params).is_err());
    }

    #[test]
    fn rejects_structure_with_no_heavy_atoms() {
        // A PDB containing only a hydrogen.
        let pdb =
            "ATOM      1  H   ALA A   1       0.000   0.000   0.000  1.00  0.00           H\nEND\n";
        let s = load_structure(pdb, "x").unwrap();
        assert!(detect_pockets(&s, &PocketParams::default()).is_err());
    }

    #[test]
    fn finds_the_cavity_in_a_hollow_box() {
        let s = load_structure(&hollow_box_pdb(), "x").unwrap();
        // The box centre is at (3*3, 3*3, 3*3) = (9, 9, 9).
        let params = PocketParams {
            spacing: 1.5,
            min_buried_directions: 10,
            probe_reach: 9.0,
            min_cluster_points: 3,
        };
        let pockets = detect_pockets(&s, &params).unwrap();
        assert!(!pockets.is_empty(), "should detect the interior cavity");
        // The best pocket's centre should be near the box centre.
        let best = &pockets[0];
        assert!(
            (best.center - Vector3::new(9.0, 9.0, 9.0)).norm() < 4.0,
            "pocket centre {:?} should be near the box centre",
            best.center
        );
        assert!(best.volume > 0.0);
        assert!((0.0..=1.0).contains(&best.score));
        assert!(!best.lining_residues.is_empty());
    }

    #[test]
    fn pocket_to_grid_box_is_dockable() {
        let s = load_structure(&hollow_box_pdb(), "x").unwrap();
        let params = PocketParams {
            spacing: 1.5,
            min_buried_directions: 10,
            probe_reach: 9.0,
            min_cluster_points: 3,
        };
        let pockets = detect_pockets(&s, &params).unwrap();
        let gb = pockets[0].to_grid_box(2.0).unwrap();
        // The grid box is centred on the pocket and has a sane edge.
        assert!(gb.size.x >= 6.0);
        assert!((gb.center - pockets[0].center).norm() < 1e-9);
    }

    #[test]
    fn solid_protein_has_few_or_no_open_cavities() {
        // A fully solid 5×5×5 block — no interior cavity to find with
        // a strict burial threshold (interior points are inside atoms).
        let mut s = String::new();
        let mut serial = 1;
        for ix in 0..5 {
            for iy in 0..5 {
                for iz in 0..5 {
                    let x = ix as f64 * 2.0;
                    let y = iy as f64 * 2.0;
                    let z = iz as f64 * 2.0;
                    s.push_str(&format!(
                        "ATOM  {serial:>5}  CA  ALA A{serial:>4}    {x:8.3}{y:8.3}{z:8.3}  1.00  0.00           C\n"
                    ));
                    serial += 1;
                }
            }
        }
        s.push_str("END\n");
        let structure = load_structure(&s, "x").unwrap();
        // With a large min_cluster_points, a solid block yields no
        // pocket (the interstitial gaps are too small).
        let params = PocketParams {
            spacing: 1.0,
            min_buried_directions: 12,
            probe_reach: 6.0,
            min_cluster_points: 50,
        };
        let pockets = detect_pockets(&structure, &params).unwrap();
        assert!(
            pockets.is_empty(),
            "solid block should expose no large cavity"
        );
    }
}
