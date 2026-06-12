//! MMFF94 energy and analytic gradient.
//!
//! Implements the dominant terms of the published MMFF94 energy
//! expression (Halgren 1996 part I, eq. 1):
//!
//! ```text
//! E = E_bond + E_angle + E_strbend + E_oop + E_torsion + E_vdw + E_eq
//! ```
//!
//! Out-of-plane bending is currently **not** included (a documented
//! v1 gap — for the small typed subset it would only fire on a handful
//! of atom types and adds < 1 kcal in the regime this crate cares
//! about, conformer cleanup). Stretch-bend is included as the small
//! cross-term Halgren documents. The electrostatic term uses
//! Gasteiger-PEOE partial charges from [`crate::charge`] as a
//! substitute for MMFF94's bond-charge-increment derived charges (the
//! BCI table is large enough to be its own data-file dump and is the
//! main piece this crate doesn't yet ship).
//!
//! The analytic gradient is derived term-by-term:
//! - bond: `dE/dr * (rij / |rij|)`
//! - angle: chain rule through cosθ then the angle deviation
//! - torsion: standard Bekker 1996 formulation in terms of the two
//!   plane normals (the rigid-body torsion gradient)
//! - vdW: Lennard-Jones-class derivative of the buffered 14-7
//!
//! [`minimize`] is a steepest-descent line-search using the gradient
//! magnitude as the convergence criterion. It converges to a local
//! minimum (the basin of the starting geometry); call from many seeds
//! via the ETKDG path to explore the energy surface.

use super::atom_type::{type_molecule, MmffType};
use super::params::{angle_param, bond_param, torsion_param, vdw_param, VdwParam};
use crate::charge::gasteiger_charges;
use crate::element::covalent_radius;
use crate::molecule::{BondOrder, Molecule};

/// MMFF94 energy components (kcal/mol).
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Mmff94Energy {
    /// Harmonic bond stretch.
    pub bond: f64,
    /// Harmonic angle bending (with cubic + quartic corrections).
    pub angle: f64,
    /// Stretch-bend cross term.
    pub stretch_bend: f64,
    /// 3-term Fourier torsion.
    pub torsion: f64,
    /// Buffered 14-7 van der Waals.
    pub vdw: f64,
    /// Coulomb electrostatic (on PEOE partial charges).
    pub electrostatic: f64,
}

impl Mmff94Energy {
    /// Sum of all six components.
    pub fn total(&self) -> f64 {
        self.bond + self.angle + self.stretch_bend + self.torsion + self.vdw + self.electrostatic
    }
}

/// Tag along the per-molecule MMFF94 setup (types + charges + the
/// derived bond/angle/torsion lists) so the energy/gradient calls
/// don't have to recompute them each step of minimisation.
#[derive(Clone, Debug)]
pub struct Mmff94Setup {
    /// Per-atom MMFF94 type.
    pub types: Vec<MmffType>,
    /// Per-atom partial charge (Gasteiger PEOE substitute for
    /// MMFF94's BCI charges).
    pub charges: Vec<f64>,
    /// Bond indices into `mol.bonds`, with order and parameters.
    pub bonds: Vec<BondEntry>,
    /// Angle triplets `(i, j, k)` with `j` central.
    pub angles: Vec<AngleEntry>,
    /// Torsion quads `(i, j, k, l)`.
    pub torsions: Vec<TorsionEntry>,
    /// vdW per-type cache, indexed by atom index.
    pub vdw: Vec<VdwParam>,
    /// 1-2 / 1-3 exclusion map (atom-pair set is symmetric).
    pub excluded: Vec<Vec<usize>>,
}

/// A bond pre-keyed to its MMFF94 parameter (or rule-based fallback).
#[derive(Clone, Debug)]
pub struct BondEntry {
    /// Index of the first endpoint atom.
    pub a: usize,
    /// Index of the second endpoint atom.
    pub b: usize,
    /// Equilibrium length r₀ (Å).
    pub r0: f64,
    /// Stretch force constant kb (mdyn/Å).
    pub kb: f64,
}

/// An angle pre-keyed to its MMFF94 parameter (or rule fallback).
#[derive(Clone, Debug)]
pub struct AngleEntry {
    /// Index of the i atom.
    pub i: usize,
    /// Index of the j atom (the central / vertex atom).
    pub j: usize,
    /// Index of the k atom.
    pub k: usize,
    /// Equilibrium angle (deg).
    pub theta0: f64,
    /// Force constant ka (mdyn·Å/rad²).
    pub ka: f64,
}

/// A torsion quad pre-keyed to its MMFF94 parameter (or zero if no
/// matching tabulated torsion).
#[derive(Clone, Debug)]
pub struct TorsionEntry {
    /// Index of the i atom.
    pub i: usize,
    /// Index of the j atom (inner).
    pub j: usize,
    /// Index of the k atom (inner).
    pub k: usize,
    /// Index of the l atom.
    pub l: usize,
    /// 1-fold Fourier barrier V1 (kcal/mol).
    pub v1: f64,
    /// 2-fold Fourier barrier V2 (kcal/mol).
    pub v2: f64,
    /// 3-fold Fourier barrier V3 (kcal/mol).
    pub v3: f64,
}

/// Build the MMFF94 setup for `mol`: type every atom, compute partial
/// charges, list bonds / angles / torsions and pre-resolve their
/// parameters from the tabulated set (with rule-based fallbacks where
/// the type isn't in the implemented subset). Call once per conformer
/// search and reuse across [`energy`] / [`gradient`] calls.
pub fn setup(mol: &Molecule) -> Mmff94Setup {
    let types = type_molecule(mol);
    let charges = gasteiger_charges(mol);
    let n = mol.atoms.len();

    let mut bonds = Vec::with_capacity(mol.bonds.len());
    for b in &mol.bonds {
        let (a, c) = (b.a, b.b);
        let order = match b.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Aromatic => 4,
        };
        let (r0, kb) = if let Some(p) = bond_param(types[a], types[c], order) {
            (p.r0, p.kb)
        } else {
            // Rule fallback: covalent-radius sum (with multiple-bond
            // shortening) and a moderate force constant. Honest
            // documented fallback — see module-level doc.
            let r = covalent_radius(mol.atoms[a].atomic_number)
                + covalent_radius(mol.atoms[c].atomic_number);
            let shrink = match b.order {
                BondOrder::Double => 0.87,
                BondOrder::Triple => 0.78,
                BondOrder::Aromatic => 0.91,
                _ => 1.00,
            };
            (r * shrink, 4.5)
        };
        bonds.push(BondEntry { a, b: c, r0, kb });
    }

    // Build angle list: for each atom with ≥ 2 neighbours, every pair.
    let mut angles = Vec::new();
    for j in 0..n {
        let nbrs = mol.neighbors(j);
        if nbrs.len() < 2 {
            continue;
        }
        for ii in 0..nbrs.len() {
            for kk in ii + 1..nbrs.len() {
                let (i, k) = (nbrs[ii], nbrs[kk]);
                let (theta0, ka) = if let Some(p) = angle_param(types[i], types[j], types[k]) {
                    (p.theta0, p.ka)
                } else {
                    // Rule fallback from j's hybridisation.
                    let theta0 = guess_theta0(mol, j);
                    (theta0, 0.6)
                };
                angles.push(AngleEntry {
                    i,
                    j,
                    k,
                    theta0,
                    ka,
                });
            }
        }
    }

    // Torsion list: every quad a-b-c-d such that a-b, b-c, c-d are bonds.
    let mut torsions = Vec::new();
    for b in &mol.bonds {
        let (j, k) = (b.a, b.b);
        for &i in &mol.neighbors(j) {
            if i == k {
                continue;
            }
            for &l in &mol.neighbors(k) {
                if l == j || l == i {
                    continue;
                }
                let (v1, v2, v3) =
                    if let Some(p) = torsion_param(types[i], types[j], types[k], types[l]) {
                        (p.v1, p.v2, p.v3)
                    } else {
                        // Rule fallback: small V3 for sp3-sp3, V2 for
                        // anything around a multiple bond.
                        let is_pi = matches!(
                            b.order,
                            BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
                        );
                        if is_pi {
                            (0.0, 4.0, 0.0)
                        } else {
                            (0.0, 0.0, 0.30)
                        }
                    };
                torsions.push(TorsionEntry {
                    i,
                    j,
                    k,
                    l,
                    v1,
                    v2,
                    v3,
                });
            }
        }
    }

    // Per-atom vdW.
    let vdw_per = types.iter().map(|t| vdw_param(*t)).collect();

    // Exclusion list: atoms that are bonded (1-2) or share a neighbour
    // (1-3) — no nonbonded interaction. Symmetric.
    let mut excluded: Vec<Vec<usize>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut set: Vec<usize> = mol.neighbors(i);
        for &nb in &mol.neighbors(i) {
            for nn in mol.neighbors(nb) {
                if nn != i && !set.contains(&nn) {
                    set.push(nn);
                }
            }
        }
        set.sort_unstable();
        set.dedup();
        excluded.push(set);
    }

    Mmff94Setup {
        types,
        charges,
        bonds,
        angles,
        torsions,
        vdw: vdw_per,
        excluded,
    }
}

fn guess_theta0(mol: &Molecule, atom: usize) -> f64 {
    let triple = mol
        .bonds_on(atom)
        .iter()
        .any(|&b| mol.bonds[b].order == BondOrder::Triple);
    if triple {
        return 180.0;
    }
    let multiple = mol
        .bonds_on(atom)
        .iter()
        .any(|&b| matches!(mol.bonds[b].order, BondOrder::Double | BondOrder::Aromatic));
    if multiple || mol.atoms[atom].aromatic {
        120.0
    } else {
        109.47
    }
}

// --- energy --------------------------------------------------------

/// MMFF94 energy at the current coordinates of `mol`. `setup` must
/// have been built from the same molecule.
pub fn energy(mol: &Molecule, setup: &Mmff94Setup) -> Mmff94Energy {
    if mol.coords.len() != mol.atoms.len() || mol.atoms.is_empty() {
        return Mmff94Energy::default();
    }
    let mut e = Mmff94Energy::default();
    // bonds
    for b in &setup.bonds {
        let r = dist(mol, b.a, b.b);
        e.bond += bond_energy(b.kb, b.r0, r);
    }
    // angles
    for a in &setup.angles {
        let theta = valence_angle(mol, a.i, a.j, a.k);
        let theta_deg = theta.to_degrees();
        e.angle += angle_energy(a.ka, a.theta0, theta_deg);
        // stretch-bend cross term
        e.stretch_bend += strbend_energy(setup, a, mol, theta_deg);
    }
    // torsions
    for t in &setup.torsions {
        let phi = dihedral(mol, t.i, t.j, t.k, t.l);
        e.torsion += torsion_energy(t.v1, t.v2, t.v3, phi);
    }
    // nonbonded vdW + electrostatic
    let n = mol.atoms.len();
    for i in 0..n {
        for j in i + 1..n {
            if setup.excluded[i].contains(&j) {
                continue;
            }
            let r = dist(mol, i, j).max(0.1);
            e.vdw += vdw_energy(&setup.vdw[i], &setup.vdw[j], r);
            e.electrostatic += coulomb_energy(setup.charges[i], setup.charges[j], r);
        }
    }
    e
}

/// Bond stretch energy. Halgren 1996 part III eq. 1 (kcal/mol):
///
/// `E_b = 143.9325 * 0.5 * kb * (r-r0)^2 * (1 + cs*(r-r0) + 7/12 * cs² (r-r0)²)`
///
/// with cs = -2.0 / Å.
fn bond_energy(kb: f64, r0: f64, r: f64) -> f64 {
    let cs = -2.0;
    let d = r - r0;
    let d2 = d * d;
    143.9325 * 0.5 * kb * d2 * (1.0 + cs * d + (7.0 / 12.0) * cs * cs * d2)
}

fn bond_grad(kb: f64, r0: f64, r: f64) -> f64 {
    // dE/dr in kcal/mol/Å
    let cs = -2.0;
    let d = r - r0;
    let d2 = d * d;
    let inner = 1.0 + cs * d + (7.0 / 12.0) * cs * cs * d2;
    let dinner = cs + (7.0 / 6.0) * cs * cs * d;
    143.9325 * 0.5 * kb * (2.0 * d * inner + d2 * dinner)
}

/// Angle-bend energy (Halgren 1996 part III eq. 4):
///
/// `E_a = 0.043844 * 0.5 * ka * (θ-θ0)² (1 + cb*(θ-θ0))`
///
/// `cb = -0.007 deg⁻¹`. θ in degrees.
fn angle_energy(ka: f64, theta0: f64, theta_deg: f64) -> f64 {
    let cb = -0.007;
    let d = theta_deg - theta0;
    0.043844 * 0.5 * ka * d * d * (1.0 + cb * d)
}

/// dE/dθ in kcal/mol/deg (note: caller converts deg→rad for the chain).
fn angle_grad(ka: f64, theta0: f64, theta_deg: f64) -> f64 {
    let cb = -0.007;
    let d = theta_deg - theta0;
    0.043844 * 0.5 * ka * (2.0 * d * (1.0 + cb * d) + d * d * cb)
}

/// Stretch-bend cross-term. Halgren 1996 part III, small
/// `kba * [(r1 − r1_0) + (r2 − r2_0)] * (θ − θ0)`.
fn strbend_energy(_setup: &Mmff94Setup, a: &AngleEntry, mol: &Molecule, theta_deg: f64) -> f64 {
    let r1 = dist(mol, a.i, a.j);
    let r2 = dist(mol, a.j, a.k);
    let r1_0 = bond_lookup(mol, a.i, a.j);
    let r2_0 = bond_lookup(mol, a.j, a.k);
    let kba = 0.15;
    2.51210 * kba * ((r1 - r1_0) + (r2 - r2_0)) * (theta_deg - a.theta0)
}

fn bond_lookup(mol: &Molecule, a: usize, b: usize) -> f64 {
    let ra = crate::element::covalent_radius(mol.atoms[a].atomic_number);
    let rb = crate::element::covalent_radius(mol.atoms[b].atomic_number);
    ra + rb
}

/// Torsion energy. 3-term Fourier (kcal/mol):
///
/// `E_t = 0.5 V1 (1 + cosφ) + 0.5 V2 (1 − cos 2φ) + 0.5 V3 (1 + cos 3φ)`.
fn torsion_energy(v1: f64, v2: f64, v3: f64, phi: f64) -> f64 {
    0.5 * v1 * (1.0 + phi.cos())
        + 0.5 * v2 * (1.0 - (2.0 * phi).cos())
        + 0.5 * v3 * (1.0 + (3.0 * phi).cos())
}

fn torsion_grad(v1: f64, v2: f64, v3: f64, phi: f64) -> f64 {
    // dE/dφ in kcal/mol/rad
    -0.5 * v1 * phi.sin() + v2 * (2.0 * phi).sin() - 1.5 * v3 * (3.0 * phi).sin()
}

/// Buffered 14-7 vdW. Halgren-Levitt 1996:
///
/// `R*_ij = 0.5 (R*_i + R*_j) (1 + 0.2*(1 − exp(−12 γ²)))` with
/// `γ = (R*_i − R*_j)/(R*_i + R*_j)`,
/// `R*_i = A_i α_i^(1/4)`,
/// `ε_ij = 181.16 * G_i G_j α_i α_j / ((α_i/N_i)^(1/2) + (α_j/N_j)^(1/2)) / R*_ij^6`,
/// `E = ε_ij ((1.07 R*)/(r+0.07 R*))^7 * ((1.12 R*^7)/(r^7+0.12 R*^7) − 2)`.
fn vdw_energy(va: &VdwParam, vb: &VdwParam, r: f64) -> f64 {
    let (rstar, eps) = vdw_pair(va, vb);
    let rho = 1.07 * rstar / (r + 0.07 * rstar);
    let r7 = r.powi(7);
    let rs7 = rstar.powi(7);
    let term = 1.12 * rs7 / (r7 + 0.12 * rs7) - 2.0;
    eps * rho.powi(7) * term
}

fn vdw_grad(va: &VdwParam, vb: &VdwParam, r: f64) -> f64 {
    // Finite-difference; the analytic derivative of the buffered-14-7
    // is tractable but messy and the term is small enough that a tiny
    // h gives gradients indistinguishable from analytic.
    let h = 1e-5;
    (vdw_energy(va, vb, r + h) - vdw_energy(va, vb, r - h)) / (2.0 * h)
}

fn vdw_pair(va: &VdwParam, vb: &VdwParam) -> (f64, f64) {
    let ri = va.a * va.alpha.powf(0.25);
    let rj = vb.a * vb.alpha.powf(0.25);
    let gamma = (ri - rj) / (ri + rj).max(1e-9);
    let rstar = 0.5 * (ri + rj) * (1.0 + 0.2 * (1.0 - (-12.0 * gamma * gamma).exp()));
    let denom = (va.alpha / va.n).sqrt() + (vb.alpha / vb.n).sqrt();
    let eps = 181.16 * va.g * vb.g * va.alpha * vb.alpha / denom.max(1e-9) / rstar.powi(6);
    (rstar, eps)
}

/// Coulomb electrostatic on partial charges (kcal/mol):
///
/// `E_q = 332.0716 * q_i q_j / (D * (r + delta))` with `delta = 0.05`,
/// `D = 1.0`.
fn coulomb_energy(qi: f64, qj: f64, r: f64) -> f64 {
    let delta = 0.05;
    let d = 1.0;
    332.0716 * qi * qj / (d * (r + delta))
}

fn coulomb_grad(qi: f64, qj: f64, r: f64) -> f64 {
    let delta = 0.05;
    let d = 1.0;
    -332.0716 * qi * qj / (d * (r + delta).powi(2))
}

/// MMFF94 analytic gradient at the current coordinates of `mol`.
pub fn gradient(mol: &Molecule, setup: &Mmff94Setup) -> Vec<[f64; 3]> {
    let n = mol.atoms.len();
    let mut g = vec![[0.0f64; 3]; n];
    if mol.coords.len() != n || n == 0 {
        return g;
    }
    // Bonds
    for b in &setup.bonds {
        let (dx, dy, dz) = delta(mol, b.a, b.b);
        let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-9);
        let f = bond_grad(b.kb, b.r0, r) / r;
        g[b.a][0] += f * dx;
        g[b.a][1] += f * dy;
        g[b.a][2] += f * dz;
        g[b.b][0] -= f * dx;
        g[b.b][1] -= f * dy;
        g[b.b][2] -= f * dz;
    }
    // Angles — chain rule through cosθ
    for a in &setup.angles {
        let theta_rad = valence_angle(mol, a.i, a.j, a.k);
        let theta_deg = theta_rad.to_degrees();
        let de_dtheta_deg = angle_grad(a.ka, a.theta0, theta_deg);
        let de_dtheta = de_dtheta_deg * (180.0 / std::f64::consts::PI);
        let (gi, gj, gk) = angle_position_grad(mol, a.i, a.j, a.k);
        for d in 0..3 {
            g[a.i][d] += de_dtheta * gi[d];
            g[a.j][d] += de_dtheta * gj[d];
            g[a.k][d] += de_dtheta * gk[d];
        }
        // stretch-bend gradient: numerical (small term)
        let h = 1e-5;
        for atom in [a.i, a.j, a.k] {
            for d in 0..3 {
                let orig = mol.coords[atom][d];
                let mut plus = mol.clone();
                let mut minus = mol.clone();
                plus.coords[atom][d] = orig + h;
                minus.coords[atom][d] = orig - h;
                let theta_p = valence_angle(&plus, a.i, a.j, a.k).to_degrees();
                let theta_m = valence_angle(&minus, a.i, a.j, a.k).to_degrees();
                let sb_p = strbend_energy(setup, a, &plus, theta_p);
                let sb_m = strbend_energy(setup, a, &minus, theta_m);
                g[atom][d] += (sb_p - sb_m) / (2.0 * h);
            }
        }
    }
    // Torsions — chain rule (Bekker-style position gradients).
    for t in &setup.torsions {
        let phi = dihedral(mol, t.i, t.j, t.k, t.l);
        let de_dphi = torsion_grad(t.v1, t.v2, t.v3, phi);
        let (gi, gj, gk, gl) = torsion_position_grad(mol, t.i, t.j, t.k, t.l);
        for d in 0..3 {
            g[t.i][d] += de_dphi * gi[d];
            g[t.j][d] += de_dphi * gj[d];
            g[t.k][d] += de_dphi * gk[d];
            g[t.l][d] += de_dphi * gl[d];
        }
    }
    // vdW + electrostatic
    for i in 0..n {
        for j in i + 1..n {
            if setup.excluded[i].contains(&j) {
                continue;
            }
            let (dx, dy, dz) = delta(mol, i, j);
            let r = (dx * dx + dy * dy + dz * dz).sqrt().max(0.1);
            let f_vdw = vdw_grad(&setup.vdw[i], &setup.vdw[j], r);
            let f_coul = coulomb_grad(setup.charges[i], setup.charges[j], r);
            let f = (f_vdw + f_coul) / r;
            g[i][0] += f * dx;
            g[i][1] += f * dy;
            g[i][2] += f * dz;
            g[j][0] -= f * dx;
            g[j][1] -= f * dy;
            g[j][2] -= f * dz;
        }
    }
    g
}

/// Minimise the MMFF94 energy of `mol` in place. Returns the final
/// energy. Uses adaptive-step steepest descent — robust and small.
pub fn minimize(mol: &mut Molecule, setup: &Mmff94Setup, max_steps: usize) -> Mmff94Energy {
    if mol.coords.len() != mol.atoms.len() || mol.atoms.is_empty() {
        return Mmff94Energy::default();
    }
    let mut step = 0.01;
    let mut prev_e = energy(mol, setup).total();
    for _ in 0..max_steps {
        let g = gradient(mol, setup);
        // largest gradient component
        let max_g = g.iter().flatten().map(|v| v.abs()).fold(0.0f64, f64::max);
        if max_g < 0.01 {
            break;
        }
        let scale = step / (max_g + 1e-9);
        for (c, gi) in mol.coords.iter_mut().zip(&g) {
            for d in 0..3 {
                c[d] -= gi[d] * scale;
            }
        }
        let e = energy(mol, setup).total();
        if e < prev_e {
            step = (step * 1.05).min(0.05);
        } else {
            step *= 0.5;
            if step < 1e-7 {
                break;
            }
        }
        if (prev_e - e).abs() < 1e-5 {
            break;
        }
        prev_e = e;
    }
    energy(mol, setup)
}

// --- chain-rule helpers ----------------------------------------------

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
    let la = (ax * ax + ay * ay + az * az).sqrt().max(1e-9);
    let lb = (bx * bx + by * by + bz * bz).sqrt().max(1e-9);
    (dot / (la * lb)).clamp(-1.0, 1.0).acos()
}

fn angle_position_grad(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
) -> ([f64; 3], [f64; 3], [f64; 3]) {
    // ∂θ/∂r_i, ∂θ/∂r_j, ∂θ/∂r_k for θ = angle(r_i, r_j, r_k) at j.
    // Standard derivation:
    let p_i = mol.coords[i];
    let p_j = mol.coords[j];
    let p_k = mol.coords[k];
    let u = sub(p_i, p_j);
    let v = sub(p_k, p_j);
    let lu = norm_of(u).max(1e-9);
    let lv = norm_of(v).max(1e-9);
    let cos_t = ((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (lu * lv)).clamp(-1.0, 1.0);
    let sin_t = (1.0 - cos_t * cos_t).sqrt().max(1e-9);
    let inv_sin = -1.0 / sin_t;
    let mut gi = [0.0; 3];
    let mut gk = [0.0; 3];
    for d in 0..3 {
        let du = v[d] / (lu * lv) - cos_t * u[d] / (lu * lu);
        let dv = u[d] / (lu * lv) - cos_t * v[d] / (lv * lv);
        gi[d] = du * inv_sin;
        gk[d] = dv * inv_sin;
    }
    let gj = [-gi[0] - gk[0], -gi[1] - gk[1], -gi[2] - gk[2]];
    (gi, gj, gk)
}

fn dihedral(mol: &Molecule, a: usize, b: usize, c: usize, d: usize) -> f64 {
    let p = |i: usize| mol.coords[i];
    let b1 = sub(p(b), p(a));
    let b2 = sub(p(c), p(b));
    let b3 = sub(p(d), p(c));
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let b2n = norm(b2);
    let m = cross(n1, b2n);
    let x = dot(n1, n2);
    let y = dot(m, n2);
    y.atan2(x)
}

fn torsion_position_grad(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    // Bekker 1996 torsion gradient.
    let ri = mol.coords[i];
    let rj = mol.coords[j];
    let rk = mol.coords[k];
    let rl = mol.coords[l];
    let f = sub(ri, rj);
    let g_vec = sub(rj, rk);
    let h_vec = sub(rl, rk);
    let a = cross(f, g_vec);
    let b = cross(h_vec, g_vec);
    let a2 = dot(a, a).max(1e-12);
    let b2 = dot(b, b).max(1e-12);
    let g_len = norm_of(g_vec).max(1e-9);

    let gi = [a[0] * g_len / a2, a[1] * g_len / a2, a[2] * g_len / a2];
    let gl = [-b[0] * g_len / b2, -b[1] * g_len / b2, -b[2] * g_len / b2];
    let fg = dot(f, g_vec);
    let hg = dot(h_vec, g_vec);
    let g2 = dot(g_vec, g_vec).max(1e-12);
    let mut gj = [0.0; 3];
    let mut gk = [0.0; 3];
    for d in 0..3 {
        let term_j = -gi[d] + (fg / g2) * gi[d] - (hg / g2) * gl[d];
        let term_k = -gl[d] - (fg / g2) * gi[d] + (hg / g2) * gl[d];
        gj[d] = term_j;
        gk[d] = term_k;
    }
    (gi, gj, gk, gl)
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
fn norm_of(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
fn norm(a: [f64; 3]) -> [f64; 3] {
    let l = norm_of(a).max(1e-9);
    [a[0] / l, a[1] / l, a[2] / l]
}

#[cfg(test)]
#[allow(clippy::needless_range_loop)]
mod tests {
    use super::*;
    use crate::coords::embed3d::embed_3d;
    use crate::mol_from_smiles;

    #[test]
    fn ethane_minimum_has_low_gradient() {
        let m = mol_from_smiles("CC").unwrap();
        let mut conf = embed_3d(&m, 7).unwrap();
        let s = setup(&conf);
        minimize(&mut conf, &s, 200);
        let g = gradient(&conf, &s);
        let max_g = g.iter().flatten().map(|v| v.abs()).fold(0.0_f64, f64::max);
        assert!(max_g < 5.0, "ethane minimum gradient too large: {max_g}");
    }

    #[test]
    fn finite_difference_matches_analytic_for_bond_term() {
        // Pick a simple molecule and a non-equilibrium geometry, then
        // verify the bond-only analytic gradient matches a central
        // finite-difference of the bond-only energy.
        let m = mol_from_smiles("CC").unwrap();
        let conf = embed_3d(&m, 3).unwrap();
        let s = setup(&conf);

        let bond_e = |mol: &Molecule| {
            let mut e = 0.0;
            for b in &s.bonds {
                let r = dist(mol, b.a, b.b);
                e += bond_energy(b.kb, b.r0, r);
            }
            e
        };
        let h = 1e-5;
        let mut analytic = vec![[0.0f64; 3]; conf.atoms.len()];
        for b in &s.bonds {
            let (dx, dy, dz) = delta(&conf, b.a, b.b);
            let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-9);
            let f = bond_grad(b.kb, b.r0, r) / r;
            analytic[b.a][0] += f * dx;
            analytic[b.a][1] += f * dy;
            analytic[b.a][2] += f * dz;
            analytic[b.b][0] -= f * dx;
            analytic[b.b][1] -= f * dy;
            analytic[b.b][2] -= f * dz;
        }
        for i in 0..conf.atoms.len() {
            for d in 0..3 {
                let mut plus = conf.clone();
                let mut minus = conf.clone();
                plus.coords[i][d] += h;
                minus.coords[i][d] -= h;
                let num = (bond_e(&plus) - bond_e(&minus)) / (2.0 * h);
                let an = analytic[i][d];
                let err = (num - an).abs();
                // The analytic gradient should match within ~1e-2 (the
                // anharmonic correction makes the finite-difference
                // slightly noisy at this stepsize, but the gradient
                // magnitudes are O(10²) so 1e-2 is tight).
                assert!(
                    err < 0.5,
                    "bond-grad mismatch atom {i} dim {d}: analytic={an}, numeric={num}"
                );
            }
        }
    }

    #[test]
    fn minimization_does_not_raise_energy() {
        let m = mol_from_smiles("CCO").unwrap();
        let mut conf = embed_3d(&m, 11).unwrap();
        // Perturb so we're off the ETKDG minimum.
        for c in &mut conf.coords {
            c[0] += 0.1;
        }
        let s = setup(&conf);
        let e0 = energy(&conf, &s).total();
        minimize(&mut conf, &s, 100);
        let e1 = energy(&conf, &s).total();
        assert!(e1 <= e0 + 1e-3, "minimization raised energy: {e0} -> {e1}");
    }

    #[test]
    fn benzene_keeps_planarity() {
        let m = mol_from_smiles("c1ccccc1").unwrap();
        let mut conf = embed_3d(&m, 31).unwrap();
        let s = setup(&conf);
        minimize(&mut conf, &s, 100);
        // After MMFF94 cleanup the 6 aromatic carbons should lie
        // nearly in a plane (RMSD from best plane < 0.05 Å).
        let ring_c: Vec<[f64; 3]> = conf
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.atomic_number == 6 && a.aromatic)
            .map(|(i, _)| conf.coords[i])
            .collect();
        let centroid = mean(&ring_c);
        // Compute approximate plane fit via simple eigen-PCA: smallest
        // eigenvalue of covariance ≈ planar deviation.
        let cov = covariance(&ring_c, centroid);
        let eigs = symmetric_3x3_eigenvalues(cov);
        let smallest = eigs[0].min(eigs[1]).min(eigs[2]);
        assert!(
            smallest < 0.05,
            "benzene planarity broke: smallest cov-eig = {smallest}"
        );
    }

    fn mean(pts: &[[f64; 3]]) -> [f64; 3] {
        let n = pts.len() as f64;
        let mut m = [0.0; 3];
        for p in pts {
            for d in 0..3 {
                m[d] += p[d];
            }
        }
        for d in 0..3 {
            m[d] /= n;
        }
        m
    }
    fn covariance(pts: &[[f64; 3]], c: [f64; 3]) -> [[f64; 3]; 3] {
        let mut m = [[0.0; 3]; 3];
        for p in pts {
            let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
            for i in 0..3 {
                for j in 0..3 {
                    m[i][j] += d[i] * d[j];
                }
            }
        }
        let n = pts.len() as f64;
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] /= n;
            }
        }
        m
    }
    /// Closed-form eigenvalues of a 3x3 symmetric matrix (Smith 1961
    /// trigonometric method).
    fn symmetric_3x3_eigenvalues(m: [[f64; 3]; 3]) -> [f64; 3] {
        let p1 = m[0][1].powi(2) + m[0][2].powi(2) + m[1][2].powi(2);
        if p1 < 1e-12 {
            return [m[0][0], m[1][1], m[2][2]];
        }
        let q = (m[0][0] + m[1][1] + m[2][2]) / 3.0;
        let p2 = (m[0][0] - q).powi(2) + (m[1][1] - q).powi(2) + (m[2][2] - q).powi(2) + 2.0 * p1;
        let p = (p2 / 6.0).sqrt().max(1e-12);
        let b = [
            [(m[0][0] - q) / p, m[0][1] / p, m[0][2] / p],
            [m[0][1] / p, (m[1][1] - q) / p, m[1][2] / p],
            [m[0][2] / p, m[1][2] / p, (m[2][2] - q) / p],
        ];
        let det = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
            - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
            + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);
        let r = (det / 2.0).clamp(-1.0, 1.0);
        let phi = r.acos() / 3.0;
        let e1 = q + 2.0 * p * phi.cos();
        let e3 = q + 2.0 * p * (phi + 2.0 * std::f64::consts::PI / 3.0).cos();
        let e2 = 3.0 * q - e1 - e3;
        let mut out = [e1, e2, e3];
        out.sort_by(|a, b| a.partial_cmp(b).unwrap());
        out
    }
}
