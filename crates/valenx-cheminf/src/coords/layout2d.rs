//! 2D depiction layout.
//!
//! [`compute_2d_coords`] assigns each atom an `(x, y)` position
//! suitable for drawing the molecule. The algorithm:
//!
//! 1. perceive the SSSR rings; lay each ring out as a regular polygon
//!    (fused rings share their common edge; the new ring is reflected
//!    to the open side);
//! 2. grow the acyclic chains outward from ring atoms (or from an
//!    arbitrary root for an acyclic molecule) placing each new atom at
//!    a standard bond length and ~120° from its parent;
//! 3. run a few rounds of a light repulsion relaxation so non-bonded
//!    atoms do not overlap.
//!
//! The result is stored in [`Molecule::coords`] with `z = 0` and
//! [`Molecule::coords_3d`] cleared. It is a *depiction* — readable, not
//! a chemically accurate geometry; use [`super::embed3d`] for that.

use crate::molecule::Molecule;
use crate::perceive::rings::{sssr, RingInfo};
use std::f64::consts::PI;

/// Standard 2D bond length (arbitrary depiction units).
const BOND_LEN: f64 = 1.5;

/// Compute a 2D depiction layout for `mol`, writing
/// [`Molecule::coords`] (z = 0). Idempotent.
pub fn compute_2d_coords(mol: &mut Molecule) {
    let n = mol.atoms.len();
    if n == 0 {
        mol.coords.clear();
        return;
    }
    let rings = sssr(mol);
    let mut pos = vec![[0.0f64, 0.0]; n];
    let mut placed = vec![false; n];

    // 1. lay out ring systems
    place_rings(mol, &rings, &mut pos, &mut placed);

    // 2. lay out acyclic atoms with a BFS chain growth
    place_chains(mol, &mut pos, &mut placed);

    // 3. any still-unplaced atom (a disconnected fragment) gets a slot
    let mut spare_x = 0.0;
    for i in 0..n {
        if !placed[i] {
            spare_x += 3.0;
            pos[i] = [spare_x, -3.0];
            placed[i] = true;
        }
    }

    // 4. light non-bonded repulsion to reduce overlaps
    relax(mol, &mut pos);

    mol.coords = pos.iter().map(|p| [p[0], p[1], 0.0]).collect();
    mol.coords_3d = false;
}

/// Place every SSSR ring as a regular polygon, fusing shared edges.
fn place_rings(mol: &Molecule, rings: &RingInfo, pos: &mut [[f64; 2]], placed: &mut [bool]) {
    let mut ring_done = vec![false; rings.rings.len()];
    // process rings in order, each fused onto the already-placed set
    for seed in 0..rings.rings.len() {
        if ring_done[seed] {
            continue;
        }
        // first ring of this system → place as a polygon at the origin
        place_polygon(&rings.rings[seed].atoms, None, pos, placed);
        ring_done[seed] = true;
        // fuse the rest of the connected ring system iteratively
        loop {
            let mut progress = false;
            for (ri, ring) in rings.rings.iter().enumerate() {
                if ring_done[ri] {
                    continue;
                }
                // does this ring share ≥2 placed atoms with the system?
                let shared: Vec<usize> = ring
                    .atoms
                    .iter()
                    .copied()
                    .filter(|&a| placed[a])
                    .collect();
                if shared.len() >= 2 {
                    place_polygon(&ring.atoms, Some(&shared), pos, placed);
                    ring_done[ri] = true;
                    progress = true;
                }
            }
            if !progress {
                break;
            }
        }
    }
    let _ = mol;
}

/// Place a ring's atoms on a regular polygon. If `shared` is given, the
/// polygon is built so it contains those already-placed atoms (a fused
/// ring); otherwise it is centred at the origin.
fn place_polygon(
    atoms: &[usize],
    shared: Option<&[usize]>,
    pos: &mut [[f64; 2]],
    placed: &mut [bool],
) {
    let k = atoms.len();
    if k < 3 {
        return;
    }
    let radius = BOND_LEN / (2.0 * (PI / k as f64).sin());

    match shared {
        None => {
            for (idx, &a) in atoms.iter().enumerate() {
                let theta = 2.0 * PI * idx as f64 / k as f64 + PI / 2.0;
                pos[a] = [radius * theta.cos(), radius * theta.sin()];
                placed[a] = true;
            }
        }
        Some(sh) if sh.len() >= 2 => {
            // place the new ring on the far side of the shared edge
            let s0 = sh[0];
            let s1 = sh[1];
            let p0 = pos[s0];
            let p1 = pos[s1];
            let mid = [(p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0];
            // edge direction and its outward normal
            let edge = [p1[0] - p0[0], p1[1] - p0[1]];
            let elen = (edge[0] * edge[0] + edge[1] * edge[1]).sqrt().max(1e-6);
            let normal = [-edge[1] / elen, edge[0] / elen];
            // centre of the new polygon, pushed along the normal away
            // from the existing ring system
            let apothem = (radius * radius - (elen / 2.0).powi(2)).max(0.0).sqrt();
            let mut center = [
                mid[0] + normal[0] * apothem,
                mid[1] + normal[1] * apothem,
            ];
            // ensure the centre is on the empty side: if any other
            // placed atom is closer to `center` than to the mirror,
            // flip
            let mirror = [
                mid[0] - normal[0] * apothem,
                mid[1] - normal[1] * apothem,
            ];
            let crowd = |c: [f64; 2]| -> f64 {
                pos.iter()
                    .enumerate()
                    .filter(|(i, _)| placed[*i])
                    .map(|(_, p)| {
                        let dx = p[0] - c[0];
                        let dy = p[1] - c[1];
                        1.0 / (dx * dx + dy * dy + 0.1)
                    })
                    .sum()
            };
            if crowd(center) > crowd(mirror) {
                center = mirror;
            }
            // angular position of s0 around `center`
            let a0 = (p0[1] - center[1]).atan2(p0[0] - center[0]);
            // determine winding from s1
            let a1 = (p1[1] - center[1]).atan2(p1[0] - center[0]);
            let step = 2.0 * PI / k as f64;
            let dir = if normalize_angle(a1 - a0) > 0.0 {
                1.0
            } else {
                -1.0
            };
            // walk the ring atoms starting at s0
            let start = atoms.iter().position(|&a| a == s0).unwrap_or(0);
            for off in 0..k {
                let a = atoms[(start + off) % k];
                if placed[a] {
                    continue;
                }
                let theta = a0 + dir * step * off as f64;
                pos[a] = [
                    center[0] + radius * theta.cos(),
                    center[1] + radius * theta.sin(),
                ];
                placed[a] = true;
            }
        }
        _ => {}
    }
}

fn normalize_angle(mut a: f64) -> f64 {
    while a > PI {
        a -= 2.0 * PI;
    }
    while a < -PI {
        a += 2.0 * PI;
    }
    a
}

/// Grow acyclic chains outward from already-placed (ring) atoms, or
/// from a root for fully-acyclic molecules.
fn place_chains(mol: &Molecule, pos: &mut [[f64; 2]], placed: &mut [bool]) {
    let n = mol.atoms.len();
    // BFS queue of placed atoms whose neighbours still need positions
    let mut queue: std::collections::VecDeque<usize> =
        (0..n).filter(|&i| placed[i]).collect();
    // if nothing placed yet (acyclic), seed atom 0 at the origin
    if queue.is_empty() && n > 0 {
        pos[0] = [0.0, 0.0];
        placed[0] = true;
        queue.push_back(0);
    }

    while let Some(u) = queue.pop_front() {
        let unplaced: Vec<usize> = mol
            .neighbors(u)
            .into_iter()
            .filter(|&v| !placed[v])
            .collect();
        if unplaced.is_empty() {
            continue;
        }
        // reference direction: away from the placed neighbour mean
        let placed_nbrs: Vec<[f64; 2]> = mol
            .neighbors(u)
            .into_iter()
            .filter(|&v| placed[v])
            .map(|v| pos[v])
            .collect();
        let base_angle = if placed_nbrs.is_empty() {
            0.0
        } else {
            let mut mx = 0.0;
            let mut my = 0.0;
            for p in &placed_nbrs {
                mx += p[0] - pos[u][0];
                my += p[1] - pos[u][1];
            }
            // point away from the mean of placed neighbours
            (my).atan2(mx) + PI
        };
        let spread = 2.0 * PI / 3.0; // ~120°
        let count = unplaced.len();
        for (idx, &v) in unplaced.iter().enumerate() {
            let offset = if count == 1 {
                0.0
            } else {
                spread * (idx as f64 / (count as f64 - 1.0) - 0.5)
            };
            let theta = base_angle + offset;
            pos[v] = [
                pos[u][0] + BOND_LEN * theta.cos(),
                pos[u][1] + BOND_LEN * theta.sin(),
            ];
            placed[v] = true;
            queue.push_back(v);
        }
    }
}

/// A few rounds of inverse-square repulsion to push overlapping atoms
/// apart, while a spring on each bond keeps lengths near `BOND_LEN`.
fn relax(mol: &Molecule, pos: &mut [[f64; 2]]) {
    let n = pos.len();
    for _ in 0..40 {
        let mut force = vec![[0.0f64, 0.0]; n];
        // repulsion between every pair
        for i in 0..n {
            for j in i + 1..n {
                let dx = pos[i][0] - pos[j][0];
                let dy = pos[i][1] - pos[j][1];
                let d2 = dx * dx + dy * dy + 0.01;
                let f = 0.4 / d2;
                let d = d2.sqrt();
                force[i][0] += f * dx / d;
                force[i][1] += f * dy / d;
                force[j][0] -= f * dx / d;
                force[j][1] -= f * dy / d;
            }
        }
        // bond springs
        for b in &mol.bonds {
            let dx = pos[b.b][0] - pos[b.a][0];
            let dy = pos[b.b][1] - pos[b.a][1];
            let d = (dx * dx + dy * dy).sqrt().max(1e-6);
            let stretch = d - BOND_LEN;
            let f = 0.3 * stretch;
            force[b.a][0] += f * dx / d;
            force[b.a][1] += f * dy / d;
            force[b.b][0] -= f * dx / d;
            force[b.b][1] -= f * dy / d;
        }
        for i in 0..n {
            // small step, capped, so the layout settles smoothly
            let step = 0.1;
            pos[i][0] += (force[i][0] * step).clamp(-0.3, 0.3);
            pos[i][1] += (force[i][1] * step).clamp(-0.3, 0.3);
        }
    }
}

/// Mean bond length of the current 2D layout — a quick quality probe.
pub fn mean_bond_length_2d(mol: &Molecule) -> f64 {
    if mol.coords.is_empty() || mol.bonds.is_empty() {
        return 0.0;
    }
    let sum: f64 = mol
        .bonds
        .iter()
        .map(|b| {
            let p = mol.coords[b.a];
            let q = mol.coords[b.b];
            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt()
        })
        .sum();
    sum / mol.bonds.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn assigns_coordinates() {
        let mut m = mol_from_smiles("CCO").unwrap();
        compute_2d_coords(&mut m);
        assert_eq!(m.coords.len(), m.atom_count());
        assert!(!m.coords_3d);
        // all z components are zero
        assert!(m.coords.iter().all(|c| c[2] == 0.0));
    }

    #[test]
    fn bonds_have_reasonable_length() {
        let mut m = mol_from_smiles("CCCCCC").unwrap();
        compute_2d_coords(&mut m);
        let mean = mean_bond_length_2d(&m);
        assert!(mean > 0.8 && mean < 2.5, "mean 2D bond length = {mean}");
    }

    #[test]
    fn ring_is_laid_out() {
        let mut m = mol_from_smiles("c1ccccc1").unwrap();
        compute_2d_coords(&mut m);
        // the 6 ring atoms should be spread out, not coincident
        let mut max_d: f64 = 0.0;
        for i in 0..6 {
            for j in i + 1..6 {
                let p = m.coords[i];
                let q = m.coords[j];
                let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt();
                max_d = max_d.max(d);
            }
        }
        // ring diameter should be ~2x the bond length
        assert!(max_d > 1.5, "benzene ring too small, max dist = {max_d}");
    }

    #[test]
    fn no_atoms_coincide() {
        let mut m = mol_from_smiles("CC(C)(C)C").unwrap(); // neopentane
        compute_2d_coords(&mut m);
        for i in 0..m.atom_count() {
            for j in i + 1..m.atom_count() {
                let p = m.coords[i];
                let q = m.coords[j];
                let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt();
                assert!(d > 0.3, "atoms {i},{j} overlap (d = {d})");
            }
        }
    }

    #[test]
    fn empty_molecule_safe() {
        let mut m = Molecule::new();
        compute_2d_coords(&mut m);
        assert!(m.coords.is_empty());
    }
}
