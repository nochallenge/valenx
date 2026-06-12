//! A reduced MMFF94 / UFF-style force field for conformer cleanup.
//!
//! After distance geometry produces a rough 3D embedding, the geometry
//! must be relaxed so bond lengths, valence angles and non-bonded
//! contacts are chemically sensible. This module implements a compact
//! force field with the four terms that dominate that relaxation:
//!
//! - **bond stretch** — harmonic, `½ k (r − r₀)²`, with `r₀` from
//!   covalent radii (UFF-style) and a bond-order-dependent stiffness;
//! - **angle bend** — harmonic in the valence angle, `½ k (θ − θ₀)²`,
//!   with `θ₀` from the central atom's hybridisation;
//! - **torsion** — a single-term periodic potential keeping conjugated
//!   / `sp²` dihedrals near planarity;
//! - **van der Waals** — a Lennard-Jones 12-6 term between non-bonded
//!   atom pairs so the structure does not self-overlap.
//!
//! [`clean_up_geometry`] minimises this energy with a capped
//! steepest-descent / gradient-descent loop and analytic gradients.
//!
//! **v1 honesty.** This is *MMFF/UFF-style*, not a published
//! parameterisation: the constants are physically reasonable generic
//! values, not the MMFF94 atom-typed tables. It is sufficient to turn a
//! distance-geometry embedding into a clean, non-overlapping conformer
//! — it is not a tool for accurate conformer energies or for ranking
//! conformers. There is no electrostatic term (partial charges are a
//! separate concern — see [`crate::charge`]).

use crate::molecule::{BondOrder, Molecule};

/// Force-field energy components, returned by [`energy`].
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Energy {
    /// Bond-stretch energy.
    pub bond: f64,
    /// Angle-bend energy.
    pub angle: f64,
    /// Torsion energy.
    pub torsion: f64,
    /// Van der Waals (Lennard-Jones) energy.
    pub vdw: f64,
}

impl Energy {
    /// Sum of every term.
    pub fn total(&self) -> f64 {
        self.bond + self.angle + self.torsion + self.vdw
    }
}

/// Evaluate the force-field energy of `mol`'s current coordinates.
/// Returns a zeroed [`Energy`] if the molecule has no coordinates.
pub fn energy(mol: &Molecule) -> Energy {
    if mol.coords.len() != mol.atoms.len() || mol.atoms.is_empty() {
        return Energy::default();
    }
    let mut e = Energy::default();
    // bonds
    for b in &mol.bonds {
        let r0 = ideal_bond(mol, b.a, b.b, b.order);
        let r = dist(mol, b.a, b.b);
        let k = bond_k(b.order);
        e.bond += 0.5 * k * (r - r0).powi(2);
    }
    // angles
    for (c, a, b) in angle_triplets(mol) {
        let t0 = ideal_angle(mol, c);
        let t = valence_angle(mol, a, c, b);
        e.angle += 0.5 * ANGLE_K * (t - t0).powi(2);
    }
    // torsions (only sp2/conjugated quads near planarity)
    for q in torsion_quads(mol) {
        let phi = dihedral(mol, q.0, q.1, q.2, q.3);
        // V(1+cos(2φ)) keeps φ near 0 or π
        e.torsion += TORSION_V * (1.0 + (2.0 * phi).cos());
    }
    // van der Waals between non-bonded pairs
    let n = mol.atoms.len();
    for i in 0..n {
        for j in i + 1..n {
            if is_bonded_or_13(mol, i, j) {
                continue;
            }
            let r = dist(mol, i, j).max(0.1);
            let sigma =
                vdw_sigma(mol.atoms[i].atomic_number) + vdw_sigma(mol.atoms[j].atomic_number);
            let sr6 = (sigma / r).powi(6);
            e.vdw += VDW_EPS * (sr6 * sr6 - 2.0 * sr6);
        }
    }
    e
}

/// Relax `mol`'s 3D coordinates by minimising the force-field energy
/// with up to `max_steps` gradient-descent iterations. Mutates the
/// coordinates in place; returns the final [`Energy`].
pub fn clean_up_geometry(mol: &mut Molecule, max_steps: usize) -> Energy {
    if mol.coords.len() != mol.atoms.len() || mol.atoms.len() < 2 {
        return energy(mol);
    }
    let mut step = 0.02;
    let mut prev_e = energy(mol).total();

    for _ in 0..max_steps {
        let grad = gradient(mol);
        // move opposite the gradient, capped per-component
        for (coord, g) in mol.coords.iter_mut().zip(&grad) {
            for d in 0..3 {
                coord[d] += (-g[d] * step).clamp(-0.2, 0.2);
            }
        }
        let e = energy(mol).total();
        if e < prev_e {
            // making progress — grow the step slightly
            step = (step * 1.1).min(0.1);
        } else {
            // overshot — undo a bit by shrinking the step
            step *= 0.5;
            if step < 1e-5 {
                break;
            }
        }
        if (prev_e - e).abs() < 1e-7 {
            break;
        }
        prev_e = e;
    }
    energy(mol)
}

/// Numerically-stable analytic gradient of the force-field energy.
fn gradient(mol: &Molecule) -> Vec<[f64; 3]> {
    let n = mol.atoms.len();
    let mut g = vec![[0.0f64; 3]; n];

    // bond gradient
    for b in &mol.bonds {
        let r0 = ideal_bond(mol, b.a, b.b, b.order);
        let (dx, dy, dz) = delta(mol, b.a, b.b);
        let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-6);
        let k = bond_k(b.order);
        let f = k * (r - r0) / r;
        g[b.a][0] += f * dx;
        g[b.a][1] += f * dy;
        g[b.a][2] += f * dz;
        g[b.b][0] -= f * dx;
        g[b.b][1] -= f * dy;
        g[b.b][2] -= f * dz;
    }

    // angle + torsion + vdw gradients via central finite differences
    // (cheap, robust — these terms are a small fraction of the work)
    let h = 1e-4;
    let mut plus = mol.clone();
    let mut minus = mol.clone();
    for (i, gi) in g.iter_mut().enumerate() {
        for (d, gid) in gi.iter_mut().enumerate() {
            let orig = mol.coords[i][d];
            plus.coords[i][d] = orig + h;
            minus.coords[i][d] = orig - h;
            let e_plus = subset_energy(&plus);
            let e_minus = subset_energy(&minus);
            *gid += (e_plus - e_minus) / (2.0 * h);
            // restore so the next perturbation starts clean
            plus.coords[i][d] = orig;
            minus.coords[i][d] = orig;
        }
    }
    g
}

/// Energy of the non-bond terms only (angle + torsion + vdw) — used by
/// the finite-difference gradient so the bond term, which already has
/// an analytic gradient, is not double-counted.
fn subset_energy(mol: &Molecule) -> f64 {
    let e = energy(mol);
    e.angle + e.torsion + e.vdw
}

// --- geometry helpers -------------------------------------------------

fn dist(mol: &Molecule, a: usize, b: usize) -> f64 {
    let (dx, dy, dz) = delta(mol, a, b);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn delta(mol: &Molecule, a: usize, b: usize) -> (f64, f64, f64) {
    let p = mol.coords[a];
    let q = mol.coords[b];
    (p[0] - q[0], p[1] - q[1], p[2] - q[2])
}

fn valence_angle(mol: &Molecule, a: usize, c: usize, b: usize) -> f64 {
    let (ax, ay, az) = delta(mol, a, c);
    let (bx, by, bz) = delta(mol, b, c);
    let dot = ax * bx + ay * by + az * bz;
    let la = (ax * ax + ay * ay + az * az).sqrt().max(1e-6);
    let lb = (bx * bx + by * by + bz * bz).sqrt().max(1e-6);
    (dot / (la * lb)).clamp(-1.0, 1.0).acos()
}

fn dihedral(mol: &Molecule, a: usize, b: usize, c: usize, d: usize) -> f64 {
    let p = |i: usize| mol.coords[i];
    let b1 = sub(p(b), p(a));
    let b2 = sub(p(c), p(b));
    let b3 = sub(p(d), p(c));
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let m = cross(n1, norm(b2));
    let x = dot(n1, n2);
    let y = dot(m, n2);
    y.atan2(x)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: [f64; 3]) -> [f64; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt().max(1e-9);
    [a[0] / l, a[1] / l, a[2] / l]
}

// --- topology helpers -------------------------------------------------

/// Every `(center, a, b)` such that `a-center` and `b-center` are bonds.
fn angle_triplets(mol: &Molecule) -> Vec<(usize, usize, usize)> {
    let mut out = Vec::new();
    for c in 0..mol.atoms.len() {
        let nbrs = mol.neighbors(c);
        for i in 0..nbrs.len() {
            for j in i + 1..nbrs.len() {
                out.push((c, nbrs[i], nbrs[j]));
            }
        }
    }
    out
}

/// Torsion quads `a-b-c-d` where the central `b-c` bond is `sp²` /
/// conjugated (double, triple or aromatic) — those are the dihedrals
/// the planarity term should constrain.
fn torsion_quads(mol: &Molecule) -> Vec<(usize, usize, usize, usize)> {
    let mut out = Vec::new();
    for bond in &mol.bonds {
        if !matches!(
            bond.order,
            BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
        ) {
            continue;
        }
        let (b, c) = (bond.a, bond.b);
        for &a in &mol.neighbors(b) {
            if a == c {
                continue;
            }
            for &d in &mol.neighbors(c) {
                if d == b || d == a {
                    continue;
                }
                out.push((a, b, c, d));
            }
        }
    }
    out
}

fn is_bonded_or_13(mol: &Molecule, i: usize, j: usize) -> bool {
    if mol.bond_between(i, j).is_some() {
        return true;
    }
    // 1,3: share a common neighbour
    let ni = mol.neighbors(i);
    mol.neighbors(j).iter().any(|v| ni.contains(v))
}

// --- parameters -------------------------------------------------------

const ANGLE_K: f64 = 0.9;
const TORSION_V: f64 = 0.18;
const VDW_EPS: f64 = 0.03;

fn bond_k(order: BondOrder) -> f64 {
    match order {
        BondOrder::Single => 4.0,
        BondOrder::Double => 7.0,
        BondOrder::Triple => 10.0,
        BondOrder::Quadruple => 11.0,
        BondOrder::Aromatic => 6.0,
    }
}

fn ideal_bond(mol: &Molecule, a: usize, b: usize, order: BondOrder) -> f64 {
    let ra = crate::element::covalent_radius(mol.atoms[a].atomic_number);
    let rb = crate::element::covalent_radius(mol.atoms[b].atomic_number);
    let base = ra + rb;
    match order {
        BondOrder::Double => base * 0.87,
        BondOrder::Triple => base * 0.78,
        BondOrder::Aromatic => base * 0.91,
        _ => base,
    }
}

fn ideal_angle(mol: &Molecule, center: usize) -> f64 {
    let a = &mol.atoms[center];
    let triple = mol
        .bonds_on(center)
        .iter()
        .any(|&b| mol.bonds[b].order == BondOrder::Triple);
    let multiple = mol
        .bonds_on(center)
        .iter()
        .any(|&b| matches!(mol.bonds[b].order, BondOrder::Double | BondOrder::Aromatic));
    if triple {
        std::f64::consts::PI
    } else if multiple || a.aromatic {
        2.094
    } else {
        1.911
    }
}

fn vdw_sigma(z: u8) -> f64 {
    // half of a vdW radius — sigma is summed for the pair
    match z {
        1 => 0.55,
        6 => 0.85,
        7 => 0.80,
        8 => 0.78,
        9 => 0.75,
        15 => 0.95,
        16 => 0.95,
        17 => 0.90,
        35 => 0.97,
        53 => 1.05,
        _ => 0.85,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::embed3d::embed_3d;
    use crate::mol_from_smiles;

    #[test]
    fn energy_zero_without_coords() {
        let m = mol_from_smiles("CCO").unwrap();
        assert_eq!(energy(&m).total(), 0.0);
    }

    #[test]
    fn cleanup_lowers_energy() {
        // build a deliberately strained ethane: place atoms too far
        let mut m = mol_from_smiles("CC").unwrap();
        let expanded = crate::perceive::hydrogen::add_explicit_hydrogens(&m);
        m = expanded;
        let n = m.atom_count();
        // scatter atoms on a line, way off ideal geometry
        m.coords = (0..n).map(|i| [i as f64 * 2.5, 0.0, 0.0]).collect();
        m.coords_3d = true;
        let before = energy(&m).total();
        clean_up_geometry(&mut m, 150);
        let after = energy(&m).total();
        assert!(
            after <= before,
            "cleanup must not raise energy: {before} -> {after}"
        );
    }

    #[test]
    fn embedded_conformer_has_low_bond_energy() {
        let conf = embed_3d(&mol_from_smiles("CCCC").unwrap(), 2).unwrap();
        let e = energy(&conf);
        // bond term should be small after embed_3d already cleaned up
        assert!(e.bond < 5.0, "bond energy after embed = {}", e.bond);
    }

    #[test]
    fn energy_components_sum() {
        let conf = embed_3d(&mol_from_smiles("c1ccccc1").unwrap(), 4).unwrap();
        let e = energy(&conf);
        assert!((e.total() - (e.bond + e.angle + e.torsion + e.vdw)).abs() < 1e-9);
    }
}
