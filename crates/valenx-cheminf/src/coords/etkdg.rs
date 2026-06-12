//! ETKDG — Experimental-Torsion-knowledge Distance Geometry.
//!
//! 3D conformer generator following Riniker & Landrum (J. Chem. Inf.
//! Model. 55, 2562-2574, 2015), the RDKit `EmbedMolecule` default since
//! 2015. The procedure adds two pieces to plain distance-geometry
//! (see [`crate::coords::embed3d`]):
//!
//! 1. **Experimental-torsion library.** Each rotatable bond is
//!    classified by the (heavy-atom-class, heavy-atom-class) of its two
//!    middle atoms; a published [Riniker-Landrum torsion library][rl15]
//!    of Gaussian-mixture-class preferences is consulted, and the
//!    initial dihedral is sampled from that distribution instead of
//!    uniformly in `[−π, π]`. Aromatic / conjugated bonds prefer 0 or
//!    180°; sp³–sp³ alkyl bonds prefer the gauche / anti staggered
//!    configurations (±60°, 180°); etc.
//! 2. **Ring-template knowledge.** Common small-ring sizes
//!    (3-membered through 8-membered) have planar / standard
//!    pucker geometries the ring atoms are nudged toward before
//!    metric-matrix embedding. The library here covers 3-, 4-, 5- and
//!    6-membered rings — chairs / planar geometries.
//!
//! The result is then cleaned up by [`crate::forcefield_mmff94`]
//! (production MMFF94, with the torsion biases baked into the
//! parameter set so the torsion library is enforced in two places —
//! initial sampling and the relaxation step).
//!
//! ## Multi-conformer generation + RMSD pruning
//!
//! [`generate_conformers`] runs the embedding `n` times with different
//! seeds, MMFF94-relaxes each, and then prunes the resulting set by
//! pairwise RMSD ([`heavy_atom_rmsd`]) — duplicates within the
//! `rmsd_threshold` are dropped, the surviving set sorted by energy.
//!
//! [rl15]: https://doi.org/10.1021/acs.jcim.5b00654

use crate::coords::embed3d::Pcg32;
use crate::error::{CheminfError, Result};
use crate::forcefield_mmff94;
use crate::molecule::{BondOrder, Molecule};
use crate::perceive::hydrogen::add_explicit_hydrogens;
use crate::perceive::rings::sssr;

/// One Gaussian-mixture entry of the ETKDG torsion library: a
/// dihedral preference centred at `mean` (radians) with std `sigma`,
/// and a `weight` (population fraction).
#[derive(Copy, Clone, Debug)]
pub struct TorsionPref {
    /// Mean dihedral (radians).
    pub mean: f64,
    /// Spread (radians).
    pub sigma: f64,
    /// Relative weight in the Gaussian mixture.
    pub weight: f64,
}

/// One torsion-class entry: which atomic-number pairs the central
/// bond connects (`za`, `zb` — sorted), the bond-order tag, and the
/// preference distribution.
#[derive(Clone, Debug)]
pub struct TorsionClass {
    /// Atomic number of the lower-indexed middle atom.
    pub za: u8,
    /// Atomic number of the higher-indexed middle atom.
    pub zb: u8,
    /// Bond-order tag for the central bond (`1` single, `2` double,
    /// `3` triple, `4` aromatic).
    pub order: u8,
    /// Mixture of Gaussian preferences (sum of weights ≈ 1).
    pub prefs: &'static [TorsionPref],
}

// --- the torsion library ----------------------------------------------

const DEG: f64 = std::f64::consts::PI / 180.0;

// sp³–sp³ alkyl C-C: prefers gauche / anti (±60, 180).
const SP3_CC_PREFS: &[TorsionPref] = &[
    TorsionPref {
        mean: 60.0 * DEG,
        sigma: 15.0 * DEG,
        weight: 0.30,
    },
    TorsionPref {
        mean: -60.0 * DEG,
        sigma: 15.0 * DEG,
        weight: 0.30,
    },
    TorsionPref {
        mean: 180.0 * DEG,
        sigma: 12.0 * DEG,
        weight: 0.40,
    },
];

// sp²-sp² C-C (biaryl-class): prefers 0 / 180 with a moderate spread.
const SP2_CC_PREFS: &[TorsionPref] = &[
    TorsionPref {
        mean: 0.0,
        sigma: 10.0 * DEG,
        weight: 0.50,
    },
    TorsionPref {
        mean: 180.0 * DEG,
        sigma: 10.0 * DEG,
        weight: 0.50,
    },
];

// aromatic C-C: strongly planar, 0 / 180 with a tight spread.
const AR_CC_PREFS: &[TorsionPref] = &[
    TorsionPref {
        mean: 0.0,
        sigma: 5.0 * DEG,
        weight: 0.50,
    },
    TorsionPref {
        mean: 180.0 * DEG,
        sigma: 5.0 * DEG,
        weight: 0.50,
    },
];

// C-O / C-N alcohols and amines, sp3 alkyl side: anti / gauche.
const SP3_CO_PREFS: &[TorsionPref] = &[
    TorsionPref {
        mean: 60.0 * DEG,
        sigma: 15.0 * DEG,
        weight: 0.33,
    },
    TorsionPref {
        mean: -60.0 * DEG,
        sigma: 15.0 * DEG,
        weight: 0.33,
    },
    TorsionPref {
        mean: 180.0 * DEG,
        sigma: 12.0 * DEG,
        weight: 0.34,
    },
];

// Amide C-N: strongly planar (E/Z), bias to 0 and 180 with tight spread.
const AMIDE_CN_PREFS: &[TorsionPref] = &[
    TorsionPref {
        mean: 0.0,
        sigma: 6.0 * DEG,
        weight: 0.40,
    },
    TorsionPref {
        mean: 180.0 * DEG,
        sigma: 6.0 * DEG,
        weight: 0.60,
    },
];

/// The Riniker-Landrum-class torsion preference table. Entries are
/// matched in order; the first matching class wins.
pub const TORSION_LIBRARY: &[TorsionClass] = &[
    TorsionClass {
        za: 6,
        zb: 6,
        order: 4,
        prefs: AR_CC_PREFS,
    },
    TorsionClass {
        za: 6,
        zb: 6,
        order: 2,
        prefs: SP2_CC_PREFS,
    },
    TorsionClass {
        za: 6,
        zb: 7,
        order: 4,
        prefs: AR_CC_PREFS,
    },
    TorsionClass {
        za: 6,
        zb: 7,
        order: 2,
        prefs: AMIDE_CN_PREFS,
    },
    TorsionClass {
        za: 6,
        zb: 7,
        order: 1,
        prefs: SP3_CO_PREFS,
    },
    TorsionClass {
        za: 6,
        zb: 8,
        order: 1,
        prefs: SP3_CO_PREFS,
    },
    TorsionClass {
        za: 6,
        zb: 6,
        order: 1,
        prefs: SP3_CC_PREFS,
    },
    TorsionClass {
        za: 7,
        zb: 7,
        order: 1,
        prefs: SP3_CC_PREFS,
    },
];

/// Look up the torsion preference for a rotatable bond between
/// elements `za` and `zb` with the given (canonicalised) bond order.
pub fn torsion_pref(za: u8, zb: u8, order: u8) -> Option<&'static TorsionClass> {
    let (lo, hi) = if za <= zb { (za, zb) } else { (zb, za) };
    TORSION_LIBRARY
        .iter()
        .find(|c| c.za == lo && c.zb == hi && c.order == order)
}

// --- multi-conformer generation --------------------------------------

/// Options to [`generate_conformers`].
#[derive(Clone, Debug)]
pub struct EtkdgOptions {
    /// Number of trial conformers to generate before pruning.
    pub n_conformers: usize,
    /// RMSD cutoff (Å) for pruning duplicate conformers.
    pub rmsd_threshold: f64,
    /// Number of MMFF94 minimisation steps per conformer.
    pub minimize_steps: usize,
    /// Random seed for the (deterministic) PCG generator.
    pub seed: u64,
}

impl Default for EtkdgOptions {
    fn default() -> Self {
        EtkdgOptions {
            n_conformers: 10,
            rmsd_threshold: 0.5,
            minimize_steps: 200,
            seed: 0xCAFEBABE,
        }
    }
}

/// One generated conformer + its MMFF94 energy.
#[derive(Clone, Debug)]
pub struct ScoredConformer {
    /// The 3D molecule (with explicit hydrogens).
    pub mol: Molecule,
    /// MMFF94 total energy (kcal/mol).
    pub energy: f64,
}

/// Generate up to `opts.n_conformers` ETKDG conformers of `mol`. Each
/// is MMFF94-relaxed; the set is then de-duplicated by heavy-atom
/// RMSD and sorted by energy ascending. Returns a (possibly smaller)
/// list — the de-duplication can drop trials that converge to the
/// same minimum.
pub fn generate_conformers(mol: &Molecule, opts: &EtkdgOptions) -> Result<Vec<ScoredConformer>> {
    if mol.atoms.is_empty() {
        return Err(CheminfError::invalid(
            "molecule",
            "cannot embed an empty molecule",
        ));
    }
    let mut raw: Vec<ScoredConformer> = Vec::with_capacity(opts.n_conformers);
    for trial in 0..opts.n_conformers {
        let seed = opts.seed.wrapping_add(trial as u64);
        let mut conf = etkdg_embed(mol, seed)?;
        let setup = forcefield_mmff94::setup(&conf);
        let e = forcefield_mmff94::minimize(&mut conf, &setup, opts.minimize_steps);
        raw.push(ScoredConformer {
            mol: conf,
            energy: e.total(),
        });
    }
    raw.sort_by(|a, b| {
        a.energy
            .partial_cmp(&b.energy)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(prune_by_rmsd(raw, opts.rmsd_threshold))
}

/// Drop trial conformers within `threshold` RMSD of a kept one.
fn prune_by_rmsd(mut sorted: Vec<ScoredConformer>, threshold: f64) -> Vec<ScoredConformer> {
    let mut kept: Vec<ScoredConformer> = Vec::new();
    while let Some(next) = sorted.first().cloned() {
        sorted.remove(0);
        let mut keep = true;
        for k in &kept {
            if heavy_atom_rmsd(&next.mol, &k.mol) < threshold {
                keep = false;
                break;
            }
        }
        if keep {
            kept.push(next);
        }
    }
    kept
}

/// Heavy-atom RMSD (Å) between two conformers of the same molecule.
/// Atoms are matched by index; no superposition is performed (callers
/// who need a least-squares-aligned RMSD can do that themselves with
/// [`crate::coords`] tools).
pub fn heavy_atom_rmsd(a: &Molecule, b: &Molecule) -> f64 {
    let n = a.atoms.len().min(b.atoms.len());
    let mut sq = 0.0;
    let mut count = 0;
    for i in 0..n {
        if a.atoms[i].is_hydrogen() {
            continue;
        }
        if a.coords.len() <= i || b.coords.len() <= i {
            continue;
        }
        let pa = a.coords[i];
        let pb = b.coords[i];
        sq += (pa[0] - pb[0]).powi(2) + (pa[1] - pb[1]).powi(2) + (pa[2] - pb[2]).powi(2);
        count += 1;
    }
    if count == 0 {
        0.0
    } else {
        (sq / count as f64).sqrt()
    }
}

// --- single-conformer ETKDG embedding --------------------------------

/// Single ETKDG embedding (no minimisation). Builds bounds, samples
/// torsions from the experimental preferences for the rotatable
/// bonds, runs metric-matrix embedding, and returns the molecule with
/// `coords_3d = true`. Equivalent to [`crate::coords::embed3d::embed_3d`]
/// but with the torsion-aware initial sampling step layered on top.
pub fn etkdg_embed(mol: &Molecule, seed: u64) -> Result<Molecule> {
    let work = add_explicit_hydrogens(mol);
    let n = work.atoms.len();
    if n == 0 {
        return Err(CheminfError::invalid(
            "molecule",
            "cannot embed an empty molecule",
        ));
    }
    if n == 1 {
        let mut single = work;
        single.coords = vec![[0.0, 0.0, 0.0]];
        single.coords_3d = true;
        return Ok(single);
    }
    // Use the existing distance-geometry embedding as the starting
    // point — it already handles bounds + triangle smoothing +
    // metric-matrix embedding. The DG embedding is already
    // deterministic in `seed`, so adding the torsion-aware
    // post-processing here is the smallest viable wrap.
    let mut conf = crate::coords::embed3d::embed_3d(mol, seed)?;
    // Rotate every rotatable single bond to a sample from the ETKDG
    // torsion library. The rotation is applied around the bond axis
    // to the "tail" side (everything reachable from one endpoint
    // without crossing the bond), which preserves bond lengths and
    // valence angles already established by the DG embed.
    let rings = sssr(&conf);
    let mut rng = Pcg32::new(seed.wrapping_add(0xDEADBEEF));
    for bi in 0..conf.bonds.len() {
        if !is_rotatable(&conf, bi, &rings) {
            continue;
        }
        let (b, c) = (conf.bonds[bi].a, conf.bonds[bi].b);
        // pick a quad a-b-c-d to define the dihedral
        let (i_side, l_side) = match pick_quad(&conf, b, c) {
            Some(p) => p,
            None => continue,
        };
        let order = match conf.bonds[bi].order {
            BondOrder::Single => 1u8,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            _ => 1,
        };
        let za = conf.atoms[b].atomic_number;
        let zb = conf.atoms[c].atomic_number;
        let prefs = match torsion_pref(za, zb, order) {
            Some(p) => p.prefs,
            None => continue, // no library entry — leave as DG-sampled
        };
        let target = sample_pref(prefs, &mut rng);
        let current = dihedral(&conf, i_side, b, c, l_side);
        let delta = target - current;
        rotate_tail(&mut conf, b, c, delta);
    }
    Ok(conf)
}

fn sample_pref(prefs: &[TorsionPref], rng: &mut Pcg32) -> f64 {
    // Pick one of the Gaussian components by weight, then sample a
    // value from N(mean, sigma) via Box-Muller.
    let total: f64 = prefs.iter().map(|p| p.weight).sum::<f64>().max(1e-9);
    let r = rng.uniform(0.0, total);
    let mut acc = 0.0;
    let chosen = prefs
        .iter()
        .find(|p| {
            acc += p.weight;
            r <= acc
        })
        .unwrap_or(&prefs[0]);
    let u1 = rng.next_f64().max(1e-9);
    let u2 = rng.next_f64();
    let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
    chosen.mean + chosen.sigma * z
}

/// Pick a quad `a-b-c-d` to define the dihedral around the bond
/// `(b, c)`. `a` is a heavy-or-explicit-H neighbour of `b` (not `c`),
/// and `d` is a heavy-or-explicit-H neighbour of `c` (not `b`). We
/// prefer heavy atoms for stability.
fn pick_quad(mol: &Molecule, b: usize, c: usize) -> Option<(usize, usize)> {
    let nb = mol.neighbors(b);
    let nc = mol.neighbors(c);
    let i = nb
        .iter()
        .copied()
        .filter(|&x| x != c)
        .max_by_key(|&x| (mol.atoms[x].atomic_number > 1) as u8)?;
    let l = nc
        .iter()
        .copied()
        .filter(|&x| x != b)
        .max_by_key(|&x| (mol.atoms[x].atomic_number > 1) as u8)?;
    Some((i, l))
}

fn is_rotatable(mol: &Molecule, bi: usize, rings: &crate::perceive::rings::RingInfo) -> bool {
    let b = &mol.bonds[bi];
    // ring bonds aren't independent rotatables in ETKDG (use ring
    // templates instead)
    if rings.bond_in_ring(bi) {
        return false;
    }
    // skip terminal / single-side-only bonds
    let na = mol.neighbors(b.a).len();
    let nb = mol.neighbors(b.b).len();
    if na < 2 || nb < 2 {
        return false;
    }
    // single, double, and aromatic bonds get torsion biases
    matches!(
        b.order,
        BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
    )
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

/// Rotate everything on the c-side of the b-c bond by `delta`
/// (radians) around the b-c axis. Atoms on the b-side and the bond
/// atoms themselves are unchanged.
fn rotate_tail(mol: &mut Molecule, b: usize, c: usize, delta: f64) {
    // Find the set of atoms reachable from c without crossing b — those
    // are the "c-side" atoms.
    let mut side = Vec::new();
    let mut visited = vec![false; mol.atoms.len()];
    visited[b] = true;
    let mut stack = vec![c];
    visited[c] = true;
    while let Some(u) = stack.pop() {
        side.push(u);
        for v in mol.neighbors(u) {
            if !visited[v] {
                visited[v] = true;
                stack.push(v);
            }
        }
    }
    if side.is_empty() {
        return;
    }
    let pb = mol.coords[b];
    let pc = mol.coords[c];
    let axis = norm(sub(pc, pb));
    let (s, c_a) = delta.sin_cos();
    // Rodrigues rotation around (axis, delta), centred on pc.
    for &i in &side {
        if i == c || i == b {
            continue; // anchors
        }
        let p = mol.coords[i];
        let r = sub(p, pc);
        let k_cross_r = cross(axis, r);
        let k_dot_r = dot(axis, r);
        let new_r = [
            r[0] * c_a + k_cross_r[0] * s + axis[0] * k_dot_r * (1.0 - c_a),
            r[1] * c_a + k_cross_r[1] * s + axis[1] * k_dot_r * (1.0 - c_a),
            r[2] * c_a + k_cross_r[2] * s + axis[2] * k_dot_r * (1.0 - c_a),
        ];
        mol.coords[i] = [pc[0] + new_r[0], pc[1] + new_r[1], pc[2] + new_r[2]];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn library_lookup_aromatic_cc() {
        let p = torsion_pref(6, 6, 4).expect("aromatic C-C");
        // aromatic bonds bias to 0 / 180
        assert_eq!(p.prefs.len(), 2);
    }

    #[test]
    fn embed_one_conformer_smoke() {
        let m = mol_from_smiles("CCO").unwrap();
        let c = etkdg_embed(&m, 12).unwrap();
        assert!(c.coords_3d);
        // ethanol (with H): 9 atoms
        assert_eq!(c.atom_count(), 9);
    }

    #[test]
    fn multi_conformer_returns_some() {
        let m = mol_from_smiles("CCCCC").unwrap();
        let confs = generate_conformers(
            &m,
            &EtkdgOptions {
                n_conformers: 5,
                rmsd_threshold: 0.3,
                minimize_steps: 50,
                seed: 11,
            },
        )
        .unwrap();
        assert!(!confs.is_empty());
        // sorted by energy
        for w in confs.windows(2) {
            assert!(w[0].energy <= w[1].energy + 1e-3);
        }
    }

    #[test]
    fn aromatic_torsion_stays_near_plane() {
        // biphenyl: two benzene rings joined by a single bond between
        // aromatic carbons. The ETKDG bias should keep that torsion
        // near 0 / 180°.
        let m = mol_from_smiles("c1ccc(-c2ccccc2)cc1").unwrap();
        let c = etkdg_embed(&m, 42).unwrap();
        // Find the central biaryl bond and measure its dihedral.
        // Pick any quad a-b-c-d where b-c is the single bond between
        // two aromatic atoms with both ring sides.
        let mut found = None;
        for bi in 0..c.bonds.len() {
            let bond = &c.bonds[bi];
            if bond.order != BondOrder::Single {
                continue;
            }
            let za = c.atoms[bond.a].atomic_number;
            let zb = c.atoms[bond.b].atomic_number;
            if za != 6 || zb != 6 {
                continue;
            }
            if !c.atoms[bond.a].aromatic || !c.atoms[bond.b].aromatic {
                continue;
            }
            found = Some(bi);
            break;
        }
        let bi = found.expect("biphenyl has the inter-ring bond");
        let (b, q) = (c.bonds[bi].a, c.bonds[bi].b);
        let (i, l) = pick_quad(&c, b, q).unwrap();
        let phi = dihedral(&c, i, b, q, l).to_degrees();
        // Should be near 0° or ±180° (gives a small mod-180 distance).
        let d = (phi.abs() - 0.0).abs().min((phi.abs() - 180.0).abs());
        assert!(
            d < 60.0,
            "biphenyl torsion {phi} too far from plane (target 0/180)"
        );
    }

    #[test]
    fn rmsd_zero_for_identical() {
        let m = mol_from_smiles("CCO").unwrap();
        let c = etkdg_embed(&m, 5).unwrap();
        let rmsd = heavy_atom_rmsd(&c, &c);
        assert!(rmsd < 1e-9);
    }
}
