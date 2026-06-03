//! Newton-Raphson 3D constraint solver with Levenberg-Marquardt damping.
//!
//! The pose vector flattens every non-fixed part's transform into 6
//! variables: `[tx, ty, tz, rx, ry, rz]`. `(tx, ty, tz)` is the
//! translation; `(rx, ry, rz)` is a **Rodrigues vector** (`||r||` =
//! rotation angle, `r/||r||` = rotation axis). Rodrigues vectors are
//! the canonical unconstrained parameterization for SO(3) — no
//! quaternion-norm constraint to satisfy.
//!
//! ## Jacobian strategy
//!
//! The translation derivatives are trivial: a Coincident mate's
//! residual is `R_b·p_b + t_b - R_a·p_a - t_a`, so `∂r/∂t_b = +I` and
//! `∂r/∂t_a = -I`. The rotation derivatives are not trivial — they
//! involve the cross product of the rotation axis with the rotated
//! body-frame anchor — and the cross-product formula is exact only at
//! the identity (angle 0); away from zero, the Jacobian of the
//! exponential map kicks in.
//!
//! **For v1 we use finite differences** for both translation and
//! rotation columns. The CPU cost is `O(6 * n_parts)` extra residual
//! evaluations per iteration, which is fine for a typical 5-part
//! assembly. We can swap to analytic-rotation entries (with the
//! exponential-map Jacobian correction) in Phase 6.5 if profiling
//! says the bottleneck moved.
//!
//! In-module tests exercise Coincident + Distance mate chains end-to-
//! end; if the Jacobian is wrong, the Newton iteration would diverge
//! or stall, so test convergence is the correctness witness.

use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};
use serde::Serialize;
use std::collections::HashMap;

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::mate::{Mate, MateKind};
use crate::part::PartTransform;

/// User-facing solver outcome.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize)]
pub enum SolverStatus {
    /// All mates satisfied within tolerance.
    Converged,
    /// Reached max iterations without converging — system may be
    /// over-constrained or inconsistent.
    MaxIterations,
}

/// Tunable solver knobs.
#[derive(Copy, Clone, Debug)]
pub struct SolverConfig {
    /// Stop when `||residual||₂ < tol`.
    pub tol: f64,
    /// Initial Levenberg-Marquardt damping.
    pub lambda_init: f64,
    /// Max Newton iterations.
    pub max_iter: usize,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            tol: 1e-9,
            lambda_init: 1e-4,
            max_iter: 100,
        }
    }
}

/// Solver outcome, including diagnostic snapshots.
#[derive(Clone, Debug, Serialize)]
pub struct SolverReport {
    /// Final status.
    pub status: SolverStatus,
    /// Iterations performed.
    pub iterations: usize,
    /// Final residual L2 norm.
    pub residual_norm: f64,
    /// Diagnostics (constraint vs DOF counts).
    pub diagnostics: SolverDiagnostics,
}

/// Constraint vs DOF accounting.
#[derive(Copy, Clone, Debug, Serialize)]
pub struct SolverDiagnostics {
    /// Total scalar residual equations (sum across un-suppressed mates).
    pub n_residuals: usize,
    /// Total variables (6 × non-fixed part count).
    pub n_variables: usize,
    /// `n_residuals - n_variables`. Negative => under-constrained;
    /// zero => exactly constrained; positive => over-constrained.
    pub dof_balance: i64,
}

/// Convert a Rodrigues vector into a [`UnitQuaternion`]. Standard
/// axis-angle parameterization with `||r||` = angle, `r/||r||` =
/// axis. Returns the identity for zero-norm input (the limit case at
/// angle 0).
fn rodrigues_to_quat(r: Vector3<f64>) -> UnitQuaternion<f64> {
    let angle = r.norm();
    if angle < 1e-12 {
        return UnitQuaternion::identity();
    }
    let axis = nalgebra::Unit::new_unchecked(r / angle);
    UnitQuaternion::from_axis_angle(&axis, angle)
}

/// Convert a [`UnitQuaternion`] back into a Rodrigues vector.
fn quat_to_rodrigues(q: &UnitQuaternion<f64>) -> Vector3<f64> {
    let (axis, angle) = q
        .axis_angle()
        .unwrap_or_else(|| (nalgebra::Unit::new_unchecked(Vector3::z()), 0.0));
    axis.into_inner() * angle
}

/// Build the flat pose vector from non-fixed parts. Layout per part
/// (6 entries): `[tx, ty, tz, rx, ry, rz]` where `(rx, ry, rz)` is
/// the Rodrigues vector for the orientation.
///
/// Returns `(vec, part_ids)` where `part_ids[i]` is the assembly
/// part id at the i-th pose-vector slot. The caller scatters results
/// back via [`apply_pose`].
pub fn pose_vector(a: &Assembly) -> (Vec<f64>, Vec<usize>) {
    let mut vec = Vec::new();
    let mut ids = Vec::new();
    for p in &a.parts {
        if p.fixed {
            continue;
        }
        let t = &p.transform.translation;
        let r = quat_to_rodrigues(&p.transform.orientation);
        vec.extend_from_slice(&[t.x, t.y, t.z, r.x, r.y, r.z]);
        ids.push(p.id);
    }
    (vec, ids)
}

/// Inverse of [`pose_vector`] — write the flat pose back to the
/// non-fixed parts' transforms.
///
/// Panics if `vec.len() != 6 * (number of non-fixed parts)` — that's
/// a programmer error, not a runtime condition.
pub fn apply_pose(a: &mut Assembly, vec: &[f64]) {
    let mut i = 0;
    for p in &mut a.parts {
        if p.fixed {
            continue;
        }
        assert!(i + 6 <= vec.len(), "pose vector too short");
        let t = Vector3::new(vec[i], vec[i + 1], vec[i + 2]);
        let r = Vector3::new(vec[i + 3], vec[i + 4], vec[i + 5]);
        p.transform = PartTransform {
            translation: t,
            orientation: rodrigues_to_quat(r),
        };
        i += 6;
    }
}

/// Compute the residual for one mate using the current transforms in
/// the assembly. Writes `mate.n_residuals()` floats into `out`.
///
/// # Errors
///
/// Returns [`AssemblyError::UnknownPart`] if either of the mate's
/// referenced part ids is missing from `a`. Returns a panic-replacing
/// length-mismatch via debug assertion when `out.len() != mate.n_residuals()` —
/// callers should size the slice via `Mate::n_residuals`.
pub fn residuals(mate: &Mate, a: &Assembly, out: &mut [f64]) -> Result<(), AssemblyError> {
    debug_assert_eq!(
        out.len(),
        mate.n_residuals(),
        "residual slice length mismatch"
    );
    if out.len() != mate.n_residuals() {
        return Err(AssemblyError::BadParameter {
            name: "out",
            reason: format!(
                "residual slice has length {} but mate expects {}",
                out.len(),
                mate.n_residuals()
            ),
        });
    }
    match &mate.kind {
        MateKind::Coincident {
            part_a,
            point_a,
            part_b,
            point_b,
        } => {
            let pa = a.get_part(*part_a)?;
            let pb = a.get_part(*part_b)?;
            let wa = pa.transform.apply_point(*point_a);
            let wb = pb.transform.apply_point(*point_b);
            out[0] = wb.x - wa.x;
            out[1] = wb.y - wa.y;
            out[2] = wb.z - wa.z;
        }
        MateKind::Distance {
            part_a,
            point_a,
            part_b,
            point_b,
            target,
        } => {
            let pa = a.get_part(*part_a)?;
            let pb = a.get_part(*part_b)?;
            let wa = pa.transform.apply_point(*point_a);
            let wb = pb.transform.apply_point(*point_b);
            out[0] = (wb - wa).norm() - *target;
        }
        MateKind::Angle {
            part_a,
            vec_a,
            part_b,
            vec_b,
            target,
        } => {
            let pa = a.get_part(*part_a)?;
            let pb = a.get_part(*part_b)?;
            let va = pa.transform.apply_vector(*vec_a);
            let vb = pb.transform.apply_vector(*vec_b);
            let cos_t = (va.dot(&vb) / (va.norm() * vb.norm())).clamp(-1.0, 1.0);
            let angle = cos_t.acos();
            out[0] = angle - *target;
        }
        MateKind::Parallel {
            part_a,
            vec_a,
            part_b,
            vec_b,
        } => {
            let pa = a.get_part(*part_a)?;
            let pb = a.get_part(*part_b)?;
            let va = pa.transform.apply_vector(*vec_a);
            let vb = pb.transform.apply_vector(*vec_b);
            // Parallel ⇔ va × vb = 0. We only need 2 of the 3
            // components (the third is linearly dependent on the
            // others); take the 2 with the largest sensitivity by
            // projecting the cross onto an orthonormal basis of
            // va's plane.
            let cross = va.cross(&vb);
            let na = va.normalize();
            // Build an arbitrary orthogonal pair to na.
            let e1 = if na.x.abs() < 0.9 {
                na.cross(&Vector3::x()).normalize()
            } else {
                na.cross(&Vector3::y()).normalize()
            };
            let e2 = na.cross(&e1);
            out[0] = cross.dot(&e1);
            out[1] = cross.dot(&e2);
        }
        MateKind::Perpendicular {
            part_a,
            vec_a,
            part_b,
            vec_b,
        } => {
            let pa = a.get_part(*part_a)?;
            let pb = a.get_part(*part_b)?;
            let va = pa.transform.apply_vector(*vec_a);
            let vb = pb.transform.apply_vector(*vec_b);
            out[0] = va.dot(&vb);
        }
        MateKind::Tangent {
            part_a,
            axis_a_origin,
            axis_a_dir,
            radius_a,
            part_b,
            axis_b_origin,
            axis_b_dir,
            radius_b,
        } => {
            let pa = a.get_part(*part_a)?;
            let pb = a.get_part(*part_b)?;
            let oa = pa.transform.apply_point(*axis_a_origin);
            let da = pa.transform.apply_vector(*axis_a_dir).normalize();
            let ob = pb.transform.apply_point(*axis_b_origin);
            // Note: axis_b_dir is intentionally unused in v1 — see the
            // Tangent semantics comment below. Once we add proper
            // skew-axis distance, db will become live.
            let _db = pb.transform.apply_vector(*axis_b_dir).normalize();
            // For two parallel axes the distance is just the
            // perpendicular component of (ob - oa) to da. We don't
            // try to handle skew axes correctly in v1 — the user is
            // expected to combine Tangent with a Parallel mate.
            let delta = ob - oa;
            let perp = delta - da * delta.dot(&da);
            let dist = perp.norm();
            out[0] = dist - (radius_a + radius_b);
        }
    }
    Ok(())
}

/// Compute the full residual vector for the assembly's un-suppressed
/// mates.
///
/// # Errors
///
/// Returns [`AssemblyError::UnknownPart`] if any mate references a
/// part id that isn't present in `a`.
pub fn assemble_residuals(a: &Assembly) -> Result<DVector<f64>, AssemblyError> {
    let n = total_residuals(a);
    let mut out = DVector::zeros(n);
    let slice = out.as_mut_slice();
    let mut row = 0;
    for m in &a.mates {
        if m.suppressed {
            continue;
        }
        let k = m.n_residuals();
        residuals(m, a, &mut slice[row..row + k])?;
        row += k;
    }
    Ok(out)
}

/// Total scalar residuals across all un-suppressed mates.
fn total_residuals(a: &Assembly) -> usize {
    a.mates
        .iter()
        .filter(|m| !m.suppressed)
        .map(|m| m.n_residuals())
        .sum()
}

/// Build a `part_id → pose-vector start column` map.
fn pose_columns(a: &Assembly) -> HashMap<usize, usize> {
    let mut map = HashMap::new();
    let mut col = 0;
    for p in &a.parts {
        if p.fixed {
            continue;
        }
        map.insert(p.id, col);
        col += 6;
    }
    map
}

/// Assemble the full Jacobian via central finite differences.
///
/// Cost: `O(n_vars * n_residuals)`. For a 5-part assembly with 30
/// vars and 15 residuals that's ~900 residual evaluations per
/// iteration — cheap.
///
/// We only perturb the columns belonging to parts that participate
/// in *some* mate (zero-row optimization is left for v2; the JᵀJ
/// solve eats the all-zero columns fine).
///
/// # Errors
///
/// Returns [`AssemblyError::UnknownPart`] if any mate references a
/// part id that isn't present in `a`.
pub fn assemble_jacobian(a: &mut Assembly, step: f64) -> Result<DMatrix<f64>, AssemblyError> {
    let (pose, _ids) = pose_vector(a);
    let n_vars = pose.len();
    let n_rows = total_residuals(a);
    let mut j = DMatrix::zeros(n_rows, n_vars);
    if n_vars == 0 || n_rows == 0 {
        return Ok(j);
    }
    let mut perturbed = pose.clone();
    for col in 0..n_vars {
        perturbed[col] = pose[col] + step;
        apply_pose(a, &perturbed);
        let r_plus = assemble_residuals(a)?;
        perturbed[col] = pose[col] - step;
        apply_pose(a, &perturbed);
        let r_minus = assemble_residuals(a)?;
        for row in 0..n_rows {
            j[(row, col)] = (r_plus[row] - r_minus[row]) / (2.0 * step);
        }
        perturbed[col] = pose[col];
    }
    // Restore the original pose.
    apply_pose(a, &pose);
    Ok(j)
}

/// Single Newton step with Levenberg-Marquardt damping. Returns the
/// step in pose-vector units, or `None` if the normal-equations
/// matrix isn't positive-definite (caller should bump `lambda`).
fn newton_step(
    a: &mut Assembly,
    lambda: f64,
    fd_step: f64,
) -> Result<Option<Vec<f64>>, AssemblyError> {
    let (pose, _) = pose_vector(a);
    let n_vars = pose.len();
    if n_vars == 0 {
        return Ok(Some(Vec::new()));
    }
    let j = assemble_jacobian(a, fd_step)?;
    let r = assemble_residuals(a)?;
    let jt = j.transpose();
    let mut jtj = &jt * &j;
    for i in 0..n_vars {
        jtj[(i, i)] += lambda;
    }
    let neg_jtr = -&jt * r;
    let Some(chol) = nalgebra::linalg::Cholesky::new(jtj) else {
        return Ok(None);
    };
    Ok(Some(chol.solve(&neg_jtr).iter().copied().collect()))
}

/// Run the Newton-Raphson 3D solver. Drives every (un-suppressed)
/// mate's residuals to zero by adjusting the pose vector of the
/// non-fixed parts.
pub fn solve(a: &mut Assembly, cfg: SolverConfig) -> Result<SolverReport, AssemblyError> {
    let diagnostics = SolverDiagnostics {
        n_residuals: total_residuals(a),
        n_variables: pose_columns(a).len() * 6,
        dof_balance: total_residuals(a) as i64 - (pose_columns(a).len() as i64 * 6),
    };
    let fd_step = 1e-7;

    let mut lambda = cfg.lambda_init;
    let mut last_norm = assemble_residuals(a)?.norm();
    for iter in 1..=cfg.max_iter {
        if last_norm < cfg.tol {
            return Ok(SolverReport {
                status: SolverStatus::Converged,
                iterations: iter - 1,
                residual_norm: last_norm,
                diagnostics,
            });
        }
        let saved_pose = pose_vector(a).0;
        let Some(step) = newton_step(a, lambda, fd_step)? else {
            lambda *= 10.0;
            if lambda > 1e12 {
                return Ok(SolverReport {
                    status: SolverStatus::MaxIterations,
                    iterations: iter,
                    residual_norm: last_norm,
                    diagnostics,
                });
            }
            continue;
        };
        let new_pose: Vec<f64> = saved_pose
            .iter()
            .zip(step.iter())
            .map(|(a, b)| a + b)
            .collect();
        apply_pose(a, &new_pose);
        let new_norm = assemble_residuals(a)?.norm();
        if new_norm < last_norm {
            last_norm = new_norm;
            lambda = (lambda * 0.5).max(1e-12);
        } else {
            // Reject: restore the previous pose and bump damping.
            apply_pose(a, &saved_pose);
            lambda *= 10.0;
            if lambda > 1e12 {
                return Ok(SolverReport {
                    status: SolverStatus::MaxIterations,
                    iterations: iter,
                    residual_norm: last_norm,
                    diagnostics,
                });
            }
        }
    }
    Ok(SolverReport {
        status: if last_norm < cfg.tol {
            SolverStatus::Converged
        } else {
            SolverStatus::MaxIterations
        },
        iterations: cfg.max_iter,
        residual_norm: last_norm,
        diagnostics,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::joint::{Joint, JointKind};
    use crate::mate::{Mate, MateKind};
    use crate::part::Part;
    use nalgebra::Vector3;

    fn unit_cube(name: &str) -> Part {
        Part::new(0, name, valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap())
    }

    #[test]
    fn pose_vector_omits_fixed_parts() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        a.add_part(p0);
        a.add_part(unit_cube("b"));
        let (pose, ids) = pose_vector(&a);
        assert_eq!(pose.len(), 6);
        assert_eq!(ids.len(), 1);
    }

    #[test]
    fn rodrigues_round_trip_at_zero_is_identity() {
        let q = rodrigues_to_quat(Vector3::zeros());
        let r = quat_to_rodrigues(&q);
        assert!(r.norm() < 1e-12);
    }

    #[test]
    fn rodrigues_round_trip_preserves_axis_and_angle() {
        use std::f64::consts::PI;
        let r_in = Vector3::new(1.0, 0.0, 0.0) * (PI / 3.0);
        let q = rodrigues_to_quat(r_in);
        let r_out = quat_to_rodrigues(&q);
        assert!((r_in - r_out).norm() < 1e-12, "got {r_out:?}");
    }

    #[test]
    fn apply_pose_round_trips_with_pose_vector() {
        let mut a = Assembly::new();
        let mut p = unit_cube("p");
        p.transform.translation = Vector3::new(1.0, 2.0, 3.0);
        p.transform.orientation = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), 0.5);
        a.add_part(p);
        let (pose, _) = pose_vector(&a);
        // Perturb then write back, then re-read — pose should round-trip.
        apply_pose(&mut a, &pose);
        let (pose2, _) = pose_vector(&a);
        assert_eq!(pose.len(), pose2.len());
        for (a, b) in pose.iter().zip(pose2.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn coincident_residual_is_world_delta() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let mut p1 = unit_cube("b");
        p1.transform.translation = Vector3::new(3.0, 4.0, 0.0);
        let id_b = a.add_part(p1);
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        ));
        let r = assemble_residuals(&a).unwrap();
        assert_eq!(r.len(), 3);
        // wb - wa = (3, 4, 0) - 0 = (3, 4, 0).
        assert!((r[0] - 3.0).abs() < 1e-12);
        assert!((r[1] - 4.0).abs() < 1e-12);
        assert!((r[2] - 0.0).abs() < 1e-12);
    }

    /// Task 14 — two unit cubes, one fixed at origin and the other
    /// initially offset. A single Coincident mate pins them together.
    /// After solve, the moving cube's vertex must coincide with the
    /// fixed cube's vertex.
    #[test]
    fn solver_converges_for_coincident_pair() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let mut p1 = unit_cube("b");
        p1.transform.translation = Vector3::new(3.0, 4.0, 0.0);
        let id_b = a.add_part(p1);
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        ));
        let report = solve(&mut a, SolverConfig::default()).unwrap();
        assert_eq!(report.status, SolverStatus::Converged, "{report:?}");
        assert!(report.residual_norm < 1e-6);
        let pb = a.get_part(id_b).unwrap();
        // After solve, b's local-frame (0,0,0) should be at world origin.
        let wb0 = pb.transform.apply_point(Vector3::zeros());
        assert!(wb0.norm() < 1e-6, "moving cube not at origin: {wb0:?}");
    }

    /// Task 15 — chain of 3 parts with Distance mates: a, b, c with
    /// dist(a,b)=5 and dist(b,c)=3. Part a is fixed at origin; b and
    /// c float.
    #[test]
    fn solver_converges_for_distance_chain() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let mut p1 = unit_cube("b");
        p1.transform.translation = Vector3::new(2.0, 0.0, 0.0);
        let id_b = a.add_part(p1);
        let mut p2 = unit_cube("c");
        p2.transform.translation = Vector3::new(7.0, 0.0, 0.0);
        let id_c = a.add_part(p2);

        a.add_mate(Mate::new(
            0,
            MateKind::Distance {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
                target: 5.0,
            },
        ));
        a.add_mate(Mate::new(
            0,
            MateKind::Distance {
                part_a: id_b,
                point_a: Vector3::zeros(),
                part_b: id_c,
                point_b: Vector3::zeros(),
                target: 3.0,
            },
        ));

        let report = solve(&mut a, SolverConfig::default()).unwrap();
        assert_eq!(report.status, SolverStatus::Converged, "{report:?}");
        assert!(report.residual_norm < 1e-6);

        let pa_origin = Vector3::zeros();
        let pb_origin = a
            .get_part(id_b)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        let pc_origin = a
            .get_part(id_c)
            .unwrap()
            .transform
            .apply_point(Vector3::zeros());
        assert!(((pb_origin - pa_origin).norm() - 5.0).abs() < 1e-6);
        assert!(((pc_origin - pb_origin).norm() - 3.0).abs() < 1e-6);
    }

    #[test]
    fn under_constrained_reports_negative_dof_balance() {
        // Single floating part, no mates → 6 vars, 0 residuals → -6.
        let mut a = Assembly::new();
        a.add_part(unit_cube("p"));
        let report = solve(&mut a, SolverConfig::default()).unwrap();
        assert_eq!(report.diagnostics.dof_balance, -6);
        assert_eq!(report.status, SolverStatus::Converged); // 0 residuals trivially converge
    }

    #[test]
    fn unused_kinematics_module_still_imports() {
        // Sanity check that the public kinematics API is reachable
        // (Phase 6A stub).
        let mut a = Assembly::new();
        a.add_part(unit_cube("p"));
        crate::kinematics::apply_all_joints(&mut a).unwrap();
        // Use Joint+JointKind to silence unused warnings until 6C
        // wires up the kinematics applier.
        let _ = Joint::new(
            0,
            JointKind::Fixed {
                part_a: 0,
                part_b: 0,
            },
        );
    }

    /// Regression: previously a mate referencing a non-existent part id
    /// caused `residuals()` to panic via `.unwrap()`. After the
    /// hardening pass, it must return [`AssemblyError::UnknownPart`].
    #[test]
    fn residuals_returns_typed_error_for_missing_part() {
        let mut a = Assembly::new();
        let mut p = unit_cube("p");
        p.fixed = true;
        let id = a.add_part(p);
        // Mate references id=999 which doesn't exist.
        let mate = Mate::new(
            0,
            MateKind::Coincident {
                part_a: id,
                point_a: Vector3::zeros(),
                part_b: 999,
                point_b: Vector3::zeros(),
            },
        );
        let mut out = vec![0.0; mate.n_residuals()];
        match residuals(&mate, &a, &mut out) {
            Err(AssemblyError::UnknownPart(999)) => {}
            other => panic!("expected UnknownPart(999), got {other:?}"),
        }
    }

    /// Regression: `assemble_residuals` also surfaces the typed error
    /// (used to panic via `residuals()`'s `.unwrap()`).
    #[test]
    fn assemble_residuals_returns_typed_error_for_missing_part() {
        let mut a = Assembly::new();
        let mut p = unit_cube("p");
        p.fixed = true;
        let id = a.add_part(p);
        a.add_mate(Mate::new(
            0,
            MateKind::Distance {
                part_a: id,
                point_a: Vector3::zeros(),
                part_b: 12345,
                point_b: Vector3::zeros(),
                target: 1.0,
            },
        ));
        match assemble_residuals(&a) {
            Err(AssemblyError::UnknownPart(12345)) => {}
            other => panic!("expected UnknownPart(12345), got {other:?}"),
        }
    }

    /// Regression: `solve()` propagates the typed error from
    /// `assemble_residuals` rather than panicking.
    #[test]
    fn solve_returns_typed_error_for_missing_part() {
        let mut a = Assembly::new();
        let mut p = unit_cube("p");
        p.fixed = true;
        let id = a.add_part(p);
        a.add_part(unit_cube("q"));
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id,
                point_a: Vector3::zeros(),
                part_b: 77777,
                point_b: Vector3::zeros(),
            },
        ));
        match solve(&mut a, SolverConfig::default()) {
            Err(AssemblyError::UnknownPart(77777)) => {}
            other => panic!("expected UnknownPart(77777), got {other:?}"),
        }
    }
}
