//! BFGS local minimization over (translation, axis-angle rotation, torsions).

use nalgebra::{UnitQuaternion, Vector3};

use crate::pose::Pose;

/// Pose serialized as a flat vector: [tx, ty, tz, rx, ry, rz, t0..tN].
/// Rotation is axis-angle (Rodrigues) so the vector is unconstrained.
pub fn pose_to_vec(pose: &Pose) -> Vec<f64> {
    let mut v = Vec::with_capacity(6 + pose.torsions.len());
    v.extend_from_slice(pose.translation.as_slice());
    let (axis, angle) = pose
        .orientation
        .axis_angle()
        .map(|(a, ang)| (a.into_inner(), ang))
        .unwrap_or((Vector3::x(), 0.0));
    v.push(axis.x * angle);
    v.push(axis.y * angle);
    v.push(axis.z * angle);
    v.extend_from_slice(&pose.torsions);
    v
}

/// Inverse of [`pose_to_vec`].
pub fn vec_to_pose(v: &[f64], n_torsions: usize) -> Pose {
    let translation = Vector3::new(v[0], v[1], v[2]);
    let rod = Vector3::new(v[3], v[4], v[5]);
    let angle = rod.norm();
    let orientation = if angle < 1e-9 {
        UnitQuaternion::identity()
    } else {
        let axis = nalgebra::Unit::new_normalize(rod / angle);
        UnitQuaternion::from_axis_angle(&axis, angle)
    };
    let torsions = v[6..6 + n_torsions].to_vec();
    Pose {
        translation,
        orientation,
        torsions,
    }
}

/// Central-difference numerical gradient of `f` at `x`.
/// Per-dimension step size: 1e-4 for translation/torsion, 1e-3 for rotation.
pub fn numerical_gradient<F: Fn(&[f64]) -> f64>(x: &[f64], n_torsions: usize, f: F) -> Vec<f64> {
    let n = 6 + n_torsions;
    let mut grad = vec![0.0; n];
    let mut xp = x.to_vec();
    for i in 0..n {
        let h = if (3..6).contains(&i) { 1e-3 } else { 1e-4 };
        let saved = xp[i];
        xp[i] = saved + h;
        let fp = f(&xp);
        xp[i] = saved - h;
        let fm = f(&xp);
        xp[i] = saved;
        grad[i] = (fp - fm) / (2.0 * h);
    }
    grad
}

use crate::eval::inter_score;
use crate::grid::GridBundle;
use crate::ligand::Ligand;

/// Run BFGS local minimization starting from `start` for at most
/// `max_iter` iterations; stop early when ||grad|| < `tol`.
///
/// Returns the optimized pose and its final score.
pub fn minimize_bfgs(
    ligand: &Ligand,
    start: &Pose,
    grids: &GridBundle,
    max_iter: usize,
    tol: f64,
) -> (Pose, f64) {
    let n_tor = ligand.n_torsions();
    let n = 6 + n_tor;
    let mut x = pose_to_vec(start);
    let f = |v: &[f64]| {
        let p = vec_to_pose(v, n_tor);
        inter_score(ligand, &p, grids)
    };
    let mut h = identity_matrix(n);
    let mut fx = f(&x);
    let mut g = numerical_gradient(&x, n_tor, f);
    for _ in 0..max_iter {
        if grad_norm(&g) < tol {
            break;
        }
        // Search direction d = -H * g
        let mut d = vec![0.0; n];
        for i in 0..n {
            let mut s = 0.0;
            for j in 0..n {
                s += h[i * n + j] * g[j];
            }
            d[i] = -s;
        }
        // Armijo backtracking line search.
        let mut step = 1.0;
        let mut x_new = x.clone();
        let mut fx_new = fx;
        for _ in 0..20 {
            for i in 0..n {
                x_new[i] = x[i] + step * d[i];
            }
            fx_new = f(&x_new);
            if fx_new < fx + 1e-4 * step * dot(&g, &d) {
                break;
            }
            step *= 0.5;
        }
        // Stop if the line search did not strictly improve fx
        // (including NaN, which is incomparable).
        if fx_new.partial_cmp(&fx) != Some(std::cmp::Ordering::Less) {
            break;
        }
        let g_new = numerical_gradient(&x_new, n_tor, f);
        // BFGS update of H (inverse Hessian).
        let s: Vec<f64> = (0..n).map(|i| x_new[i] - x[i]).collect();
        let y: Vec<f64> = (0..n).map(|i| g_new[i] - g[i]).collect();
        let sy = dot(&s, &y);
        if sy.abs() > 1e-10 {
            let rho = 1.0 / sy;
            let mut hy = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    hy[i] += h[i * n + j] * y[j];
                }
            }
            let yhy = dot(&y, &hy);
            let mut h_new = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    h_new[i * n + j] = h[i * n + j] + rho * s[i] * s[j]
                        - rho * (hy[i] * s[j] + s[i] * hy[j])
                        + rho * rho * yhy * s[i] * s[j];
                }
            }
            h = h_new;
        }
        x = x_new;
        fx = fx_new;
        g = g_new;
    }
    (vec_to_pose(&x, n_tor), fx)
}

fn identity_matrix(n: usize) -> Vec<f64> {
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        m[i * n + i] = 1.0;
    }
    m
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn grad_norm(g: &[f64]) -> f64 {
    g.iter().map(|x| x * x).sum::<f64>().sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_pose_roundtrips() {
        let p = Pose::identity(3);
        let v = pose_to_vec(&p);
        let p2 = vec_to_pose(&v, 3);
        assert_eq!(p.translation, p2.translation);
        // Identity quaternion may roundtrip as a slightly different
        // but equivalent representation; check rotation magnitude is zero.
        let (_, ang) = p2
            .orientation
            .axis_angle()
            .unwrap_or((Vector3::x_axis(), 0.0));
        assert!(ang.abs() < 1e-9);
        assert_eq!(p.torsions, p2.torsions);
    }

    #[test]
    fn translation_roundtrips() {
        let mut p = Pose::identity(2);
        p.translation = Vector3::new(1.5, -2.0, 3.25);
        p.torsions = vec![0.1, -0.2];
        let v = pose_to_vec(&p);
        let p2 = vec_to_pose(&v, 2);
        assert!((p.translation - p2.translation).norm() < 1e-12);
        assert!((p.torsions[0] - p2.torsions[0]).abs() < 1e-12);
        assert!((p.torsions[1] - p2.torsions[1]).abs() < 1e-12);
    }

    use crate::atom_type::Ad4AtomType;
    use crate::eval::inter_score;
    use crate::grid::GridBundle;
    use crate::ligand::Ligand;
    use crate::receptor::{Receptor, ReceptorAtom};

    fn one_carbon_each() -> (Ligand, Receptor, GridBundle) {
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        };
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C 
ENDROOT
TORSDOF 0
";
        let lig = Ligand::from_pdbqt(pdbqt).unwrap();
        let grids = GridBundle::build(
            &receptor,
            &lig,
            Vector3::new(-5.0, -5.0, -5.0),
            0.5,
            (21, 21, 21),
        );
        (lig, receptor, grids)
    }

    #[test]
    fn gradient_points_toward_attractive_well() {
        // Single C atom ligand, single C receptor at origin. Place
        // ligand at (3, 0, 0) — gradient in x should be NEGATIVE
        // (energy decreases as we move toward origin).
        let (lig, _r, grids) = one_carbon_each();
        let mut pose = Pose::identity(0);
        pose.translation = Vector3::new(3.0, 0.0, 0.0);
        let v = pose_to_vec(&pose);
        let g = numerical_gradient(&v, lig.n_torsions(), |w| {
            let p = vec_to_pose(w, lig.n_torsions());
            inter_score(&lig, &p, &grids)
        });
        assert!(g[0] < 0.0, "expected dE/dx < 0, got {}", g[0]);
    }

    #[test]
    fn bfgs_drives_pose_into_attractive_well() {
        let (lig, _r, grids) = one_carbon_each();
        let mut pose = Pose::identity(0);
        pose.translation = Vector3::new(3.0, 0.0, 0.0);
        let initial_score = inter_score(&lig, &pose, &grids);
        let (final_pose, final_score) = minimize_bfgs(&lig, &pose, &grids, 50, 1e-3);
        assert!(
            final_score < initial_score,
            "BFGS should have improved score: initial={initial_score}, final={final_score}"
        );
        // The minimum for a C-C pair is at surface distance d = 0 — i.e.
        // centre-to-centre distance equal to the sum of VDW radii
        // (1.9 + 1.9 = 3.8 Å). Starting at 3 Å, BFGS should pull the
        // ligand closer to that equilibrium separation.
        let target = 3.8;
        let initial_dist = 3.0;
        let final_dist = final_pose.translation.norm();
        assert!(
            (final_dist - target).abs() < (initial_dist - target).abs(),
            "translation should be closer to equilibrium {target} Å than initial {initial_dist} Å, got {final_dist}"
        );
    }
}
