//! A compact, self-contained classical force field for the QM/MM
//! environment: Lennard-Jones + Coulomb nonbonded interactions, in
//! atomic units (positions bohr, charge e, σ bohr, ε hartree → energy
//! hartree, force hartree/bohr).
//!
//! This is the **mechanical-embedding** coupling: the QM region is solved
//! by qchem in isolation, and its interaction with the explicit classical
//! environment is purely classical LJ + Coulomb. [`classical_forces`]
//! sums the **MM–MM** pairs (solvent–solvent) and the **QM–MM** pairs
//! (the coupling), but never the **QM–QM** pairs — those are qchem's job.
//!
//! Honest scope: a lean v1. Monatomic explicit solvent (no MM bonded
//! terms yet); reusing `valenx-md`'s full OPLS-AA force field (bonded
//! water, validated parameters) and electrostatic embedding (MM charges
//! polarizing the QM density) are the documented upgrades.

/// A classical particle, atomic units: position bohr, charge e, LJ σ
/// bohr, LJ ε hartree. Used for both QM atoms (LJ/charge assigned for the
/// coupling) and MM environment atoms.
#[derive(Clone, Debug, PartialEq)]
pub struct Particle {
    /// Position in bohr.
    pub pos_bohr: [f64; 3],
    /// Partial charge in e (= atomic units).
    pub charge: f64,
    /// Lennard-Jones σ in bohr.
    pub sigma_bohr: f64,
    /// Lennard-Jones ε in hartree.
    pub epsilon_hartree: f64,
}

#[inline]
fn add_scaled(t: &mut [f64; 3], v: &[f64; 3], s: f64) {
    t[0] += s * v[0];
    t[1] += s * v[1];
    t[2] += s * v[2];
}

/// Lennard-Jones (Lorentz-Berthelot combining) + Coulomb energy of a
/// pair and the **force on `a`** (the force on `b` is its negation).
/// Atomic units. Coincident atoms contribute nothing (guarded).
pub fn pair_lj_coulomb(a: &Particle, b: &Particle) -> (f64, [f64; 3]) {
    let r_vec = [
        b.pos_bohr[0] - a.pos_bohr[0],
        b.pos_bohr[1] - a.pos_bohr[1],
        b.pos_bohr[2] - a.pos_bohr[2],
    ];
    let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
    let r = r2.sqrt();
    if r < 1e-9 {
        return (0.0, [0.0; 3]);
    }
    let inv_r = 1.0 / r;
    // r̂ = (b - a) / r.
    let rhat = [r_vec[0] * inv_r, r_vec[1] * inv_r, r_vec[2] * inv_r];

    // Coulomb: V = qa·qb / r; F_a = −qa·qb / r² · r̂.
    let e_coul = a.charge * b.charge * inv_r;
    let coeff_coul = -a.charge * b.charge * inv_r * inv_r;

    // Lennard-Jones: V = 4ε[(σ/r)¹² − (σ/r)⁶]; F_a = −24ε/r[2(σ/r)¹² − (σ/r)⁶]·r̂.
    let sigma = 0.5 * (a.sigma_bohr + b.sigma_bohr);
    let eps = (a.epsilon_hartree * b.epsilon_hartree).sqrt();
    let (e_lj, coeff_lj) = if eps > 0.0 && sigma > 0.0 {
        let sr = sigma * inv_r;
        let sr6 = sr.powi(6);
        let sr12 = sr6 * sr6;
        (
            4.0 * eps * (sr12 - sr6),
            -24.0 * eps * inv_r * (2.0 * sr12 - sr6),
        )
    } else {
        (0.0, 0.0)
    };

    let coeff = coeff_coul + coeff_lj;
    let force_on_a = [coeff * rhat[0], coeff * rhat[1], coeff * rhat[2]];
    (e_coul + e_lj, force_on_a)
}

/// Classical LJ + Coulomb forces for QM/MM **mechanical embedding**:
/// the MM–MM pairs (solvent–solvent) and the QM–MM pairs (coupling), in
/// atomic units. QM–QM pairs are excluded — qchem owns the QM region.
///
/// Returns `(energy_hartree, qm_forces, mm_forces)`, each force list
/// indexed like its input slice.
pub fn classical_forces(
    qm: &[Particle],
    mm: &[Particle],
) -> (f64, Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut energy = 0.0;
    let mut f_qm = vec![[0.0; 3]; qm.len()];
    let mut f_mm = vec![[0.0; 3]; mm.len()];

    // MM–MM (solvent–solvent).
    for i in 0..mm.len() {
        for j in (i + 1)..mm.len() {
            let (e, fa) = pair_lj_coulomb(&mm[i], &mm[j]);
            energy += e;
            add_scaled(&mut f_mm[i], &fa, 1.0);
            add_scaled(&mut f_mm[j], &fa, -1.0);
        }
    }
    // QM–MM (the coupling).
    for (qi, qp) in qm.iter().enumerate() {
        for (mi, mp) in mm.iter().enumerate() {
            let (e, fa) = pair_lj_coulomb(qp, mp);
            energy += e;
            add_scaled(&mut f_qm[qi], &fa, 1.0);
            add_scaled(&mut f_mm[mi], &fa, -1.0);
        }
    }
    (energy, f_qm, f_mm)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn charge_only(pos: [f64; 3], q: f64) -> Particle {
        Particle {
            pos_bohr: pos,
            charge: q,
            sigma_bohr: 0.0,
            epsilon_hartree: 0.0,
        }
    }

    #[test]
    fn coulomb_like_charges_repel() {
        // Two +1 charges 2 bohr apart on z: V = 1/2; force on a points −z
        // (away from b at +z), magnitude qa·qb/r² = 1/4.
        let a = charge_only([0.0, 0.0, 0.0], 1.0);
        let b = charge_only([0.0, 0.0, 2.0], 1.0);
        let (e, fa) = pair_lj_coulomb(&a, &b);
        assert!((e - 0.5).abs() < 1e-12, "energy {e}");
        assert!((fa[2] - (-0.25)).abs() < 1e-12, "force_z {}", fa[2]);
    }

    #[test]
    fn lennard_jones_attracts_beyond_the_minimum() {
        // σ = 2 bohr → r_min = 2^(1/6)·σ ≈ 2.245. At r = 4 (> r_min) the
        // force on a is attractive (toward b at +z → +z).
        let p = |z: f64| Particle {
            pos_bohr: [0.0, 0.0, z],
            charge: 0.0,
            sigma_bohr: 2.0,
            epsilon_hartree: 0.01,
        };
        let (_, fa) = pair_lj_coulomb(&p(0.0), &p(4.0));
        assert!(fa[2] > 0.0, "expected attraction (+z), got {}", fa[2]);
    }

    #[test]
    fn newtons_third_law_holds() {
        let a = Particle { pos_bohr: [0.0, 0.0, 0.0], charge: 0.5, sigma_bohr: 2.0, epsilon_hartree: 0.01 };
        let b = Particle { pos_bohr: [0.3, 0.7, 1.6], charge: -0.4, sigma_bohr: 2.5, epsilon_hartree: 0.02 };
        let (_, fa) = pair_lj_coulomb(&a, &b);
        let (_, fb) = pair_lj_coulomb(&b, &a);
        for d in 0..3 {
            assert!((fa[d] + fb[d]).abs() < 1e-12, "third law dim {d}: {} vs {}", fa[d], fb[d]);
        }
    }

    #[test]
    fn classical_forces_excludes_qm_qm_pairs() {
        // Two QM particles, no MM → no classical energy (QM-QM is qchem's).
        let qm = vec![charge_only([0.0; 3], 1.0), charge_only([0.0, 0.0, 2.0], 1.0)];
        let (e, fqm, fmm) = classical_forces(&qm, &[]);
        assert_eq!(e, 0.0);
        assert!(fqm.iter().all(|f| f == &[0.0; 3]));
        assert!(fmm.is_empty());
    }

    #[test]
    fn classical_forces_couples_qm_to_mm_third_law() {
        // One QM + one MM charge → equal-and-opposite forces, energy = qq/r.
        let qm = vec![charge_only([0.0, 0.0, 0.0], 1.0)];
        let mm = vec![charge_only([0.0, 0.0, 2.0], -1.0)];
        let (e, fqm, fmm) = classical_forces(&qm, &mm);
        assert!((e - (-0.5)).abs() < 1e-12, "energy {e}");
        for d in 0..3 {
            assert!((fqm[0][d] + fmm[0][d]).abs() < 1e-12);
        }
        // Opposite charges attract: QM atom pulled +z toward the MM atom.
        assert!(fqm[0][2] > 0.0);
    }
}
