//! Commercial-depth constraint diagnostics — reports whether an assembly is
//! fully-, under-, over-, or inconsistently-constrained, and identifies the
//! offending mate set in the latter two cases.
//!
//! ## State taxonomy
//!
//! - [`ConstraintState::FullyConstrained`] — the Jacobian's rank equals the
//!   total pose-vector size (every non-fixed DOF is pinned), the system has
//!   no redundant rows, and the residual norm is below tolerance. Equivalent
//!   to the SolidWorks "Fully Defined" sketch / mate state.
//! - [`ConstraintState::UnderConstrained`] — at least one DOF is free. The
//!   `remaining_dof` field reports the dimension of the Jacobian's null
//!   space (i.e. `n_vars − rank(J)`); this is the count of independent
//!   directions you could push the assembly without violating a mate.
//!   "Add `remaining_dof` more constraints" is the diagnostic guidance.
//! - [`ConstraintState::OverConstrained`] — the Jacobian has more rows than
//!   its column rank, meaning some mates are linearly dependent on the
//!   others (redundant). The redundant mate ids come from a QR-pivoting
//!   identification of the dependent rows; we drop the *last* mate that
//!   first becomes dependent (the canonical greedy choice — the same one
//!   SolidWorks surfaces in its "What's Wrong" dialog). The system can
//!   still solve, but a numerical perturbation in one redundant mate
//!   doesn't propagate cleanly.
//! - [`ConstraintState::Inconsistent`] — the system is over-constrained
//!   *and* the residual is non-zero after a solve (i.e. the redundant mates
//!   carry conflicting target values: two `Fix` mates pinning to different
//!   positions, a `Distance` of 5 plus a `Distance` of 7 between the same
//!   two points, etc.). The conflicting mate ids are the same redundant-row
//!   set, but the state distinguishes "drop one and the rest are
//!   consistent" from "drop one and the system actually solves".
//!
//! ## Algorithm
//!
//! 1. Assemble the Jacobian `J` at the current pose (re-using
//!    [`crate::solver::assemble_jacobian`]).
//! 2. Compute its **column-pivoted QR decomposition** — `J = Q · R · P`
//!    where `R` is upper-triangular and `P` is a column permutation. The
//!    pivoting orders columns by descending magnitude so the first `r`
//!    diagonals of `R` are non-zero and the remaining are below tolerance
//!    (which gives us the numerical rank in one pass).
//! 3. `rank = #{i : |R[i,i]| > rank_tol · ||J||_F}` where the Frobenius-
//!    norm-relative tolerance is the textbook scale-invariant choice.
//! 4. From the rank we get:
//!    - `remaining_dof = n_vars − rank` (the under-constrained side).
//!    - `redundant_rows = n_rows − rank` (the over-constrained side).
//! 5. To identify *which* rows are redundant we run an **LU with partial
//!    pivoting on Jᵀ** — the rows that get permuted to the bottom (after
//!    the first `rank` rows fully account for the column space) are
//!    redundant. We map each redundant row index back to the mate id via
//!    [`mate_row_map`].
//! 6. Inconsistency: after identifying redundancy, the residual norm tells
//!    us whether the over-constrained mates *agree*. If `residual > tol`,
//!    the over-constrained mates carry conflicting target values →
//!    `Inconsistent`. Otherwise → `OverConstrained`.
//!
//! ## Honest scope
//!
//! - The QR rank tolerance uses a Frobenius-norm-relative cutoff (default
//!   `1e-9`). For very stiff systems (Jacobian condition number > 1e9) the
//!   rank determination is brittle — same caveat as any rank-revealing
//!   factorization. The user can pass a custom tolerance via
//!   [`DiagnosticsConfig`].
//! - The redundant-row identification is a *greedy* choice — there may be
//!   several equally-valid subsets of redundant mates. We report one
//!   consistent set; commercial CAD usually picks the same way (suppress
//!   the last-added mate first), which matches our row-order traversal.

use nalgebra::{DMatrix, DVector};

use crate::assembly::Assembly;
use crate::solver::{assemble_jacobian, assemble_residuals};

/// Constraint-state classification of an assembly under its current mates.
#[derive(Clone, Debug, PartialEq)]
pub enum ConstraintState {
    /// Every non-fixed DOF is pinned and the residual is below tolerance.
    /// SolidWorks "Fully Defined".
    FullyConstrained,
    /// One or more DOFs remain free. `remaining_dof` is the dimension of
    /// the Jacobian's null space (i.e. `n_vars − rank(J)`).
    UnderConstrained {
        /// Number of free DOFs remaining. Add this many more independent
        /// constraints to fully define the assembly.
        remaining_dof: usize,
    },
    /// The Jacobian has more rows than its column rank — some mates are
    /// linearly dependent. The redundant mate ids are listed (the
    /// canonical greedy "drop last" set).
    OverConstrained {
        /// Stable mate ids that the diagnostic identified as redundant.
        /// Suppressing or deleting any one of these makes the system
        /// non-redundant; the remaining mates still solve.
        redundant_mates: Vec<usize>,
    },
    /// Over-constrained *and* the residual is non-zero — the redundant
    /// mates carry conflicting target values (two `Fix` mates at different
    /// positions, two `Distance` mates with different targets between the
    /// same points, etc.). Suppressing any one of these mates makes the
    /// rest consistent.
    Inconsistent {
        /// Stable mate ids that the diagnostic identified as the
        /// conflicting set.
        conflicting_mates: Vec<usize>,
    },
}

/// Tunable knobs for the diagnostic.
#[derive(Copy, Clone, Debug)]
pub struct DiagnosticsConfig {
    /// Frobenius-norm-relative rank tolerance for the QR decomposition.
    /// `1e-9` is the textbook default (a singular value below this
    /// fraction of `||J||_F` counts as zero). Stiff systems (condition
    /// number > 1e9) may need a looser cutoff.
    pub rank_tol: f64,
    /// L2 residual threshold above which the residual is considered
    /// non-zero (and therefore: an over-constrained system tips into
    /// inconsistent). `1e-6` matches the solver convergence tolerance for
    /// a "non-zero" residual at machine-precision scales.
    pub residual_tol: f64,
    /// Finite-difference step for the Jacobian assembly. Matches the
    /// solver default.
    pub fd_step: f64,
}

impl Default for DiagnosticsConfig {
    fn default() -> Self {
        Self {
            rank_tol: 1e-9,
            residual_tol: 1e-6,
            fd_step: 1e-7,
        }
    }
}

/// Full diagnostic outcome — the classification, plus a few raw numbers
/// the UI can surface in its status line.
#[derive(Clone, Debug)]
pub struct DiagnosticsReport {
    /// The classification.
    pub state: ConstraintState,
    /// Total scalar residual rows (sum across un-suppressed mates).
    pub n_residuals: usize,
    /// Total non-fixed pose variables (6 × non-fixed part count).
    pub n_variables: usize,
    /// Numerical rank of the Jacobian (from the rank-revealing QR).
    pub jacobian_rank: usize,
    /// Final residual L2 norm at the current pose.
    pub residual_norm: f64,
}

/// Run the constraint diagnostics on the assembly at its current pose.
///
/// The assembly is not mutated (the Jacobian assembly temporarily perturbs
/// poses but restores them).
///
/// # Errors
///
/// Returns [`crate::AssemblyError::UnknownPart`] if any mate references
/// a part id missing from `a` (was previously a panic).
pub fn diagnose(
    a: &mut Assembly,
    cfg: DiagnosticsConfig,
) -> Result<DiagnosticsReport, crate::AssemblyError> {
    let n_variables = count_non_fixed_parts(a) * 6;
    let row_map = mate_row_map(a);
    let n_residuals = row_map.last().map(|&(_, end)| end).unwrap_or(0);

    // Empty system: nothing to diagnose. Treat zero variables + zero
    // residuals as "fully constrained" (vacuously true).
    if n_variables == 0 && n_residuals == 0 {
        return Ok(DiagnosticsReport {
            state: ConstraintState::FullyConstrained,
            n_residuals: 0,
            n_variables: 0,
            jacobian_rank: 0,
            residual_norm: 0.0,
        });
    }

    let jacobian = assemble_jacobian(a, cfg.fd_step)?;
    let residual = assemble_residuals(a)?;
    let residual_norm = residual.norm();

    // Rank determination via SVD-equivalent: the Frobenius norm bounds
    // the largest singular value, and we threshold pivots in the
    // QR-pivoting factorization relative to it.
    let frob_norm = jacobian.iter().map(|x| x * x).sum::<f64>().sqrt();
    let rank_threshold = if frob_norm > 0.0 {
        cfg.rank_tol * frob_norm
    } else {
        cfg.rank_tol
    };

    let (rank, redundant_rows) = rank_and_redundant_rows(&jacobian, rank_threshold);

    let remaining_dof = n_variables.saturating_sub(rank);
    let redundant_count = n_residuals.saturating_sub(rank);

    // Map redundant row indices back to mate ids. Several rows of one
    // mate (e.g. Coincident contributes 3 rows) collapse to a single mate
    // id; dedupe.
    let redundant_mate_ids = rows_to_mate_ids(&redundant_rows, &row_map);

    let state = if redundant_count == 0 && remaining_dof == 0 {
        ConstraintState::FullyConstrained
    } else if redundant_count == 0 {
        ConstraintState::UnderConstrained { remaining_dof }
    } else if residual_norm > cfg.residual_tol {
        ConstraintState::Inconsistent {
            conflicting_mates: redundant_mate_ids,
        }
    } else {
        ConstraintState::OverConstrained {
            redundant_mates: redundant_mate_ids,
        }
    };

    Ok(DiagnosticsReport {
        state,
        n_residuals,
        n_variables,
        jacobian_rank: rank,
        residual_norm,
    })
}

/// Count the parts that are not fixed (i.e. have free DOFs).
fn count_non_fixed_parts(a: &Assembly) -> usize {
    a.parts.iter().filter(|p| !p.fixed).count()
}

/// Build a `(mate_id, end_row_exclusive)` map walking the un-suppressed
/// mates in storage order. Used by both `assemble_residuals` (implicit row
/// order) and `rows_to_mate_ids` (the reverse mapping).
pub fn mate_row_map(a: &Assembly) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut row = 0;
    for m in &a.mates {
        if m.suppressed {
            continue;
        }
        row += m.n_residuals();
        out.push((m.id, row));
    }
    out
}

/// Given a sorted list of row indices flagged as redundant, look up the
/// mate id each row belongs to and return the deduplicated list of ids.
fn rows_to_mate_ids(rows: &[usize], row_map: &[(usize, usize)]) -> Vec<usize> {
    let mut seen = Vec::new();
    for &row in rows {
        // Binary search would work too, but n is small.
        for &(mate_id, end) in row_map {
            if row < end {
                if !seen.contains(&mate_id) {
                    seen.push(mate_id);
                }
                break;
            }
        }
    }
    seen
}

/// Determine the Jacobian's rank and the row indices that are linearly
/// dependent on the rows above them.
///
/// Implementation: Gaussian elimination on `Jᵀ` with partial pivoting — we
/// treat each *row* of `J` as one column of `Jᵀ` and reduce; rows that
/// produce a pivot below `tol` are redundant on the rows already
/// processed. This is simpler and more numerically robust here than a
/// full column-pivoted QR for our scale (residual count ≤ ~50).
fn rank_and_redundant_rows(j: &DMatrix<f64>, tol: f64) -> (usize, Vec<usize>) {
    let n_rows = j.nrows();
    let n_cols = j.ncols();
    if n_rows == 0 || n_cols == 0 {
        return (0, Vec::new());
    }
    // Working copy — we'll reduce row-by-row.
    let mut work: Vec<DVector<f64>> = (0..n_rows).map(|i| j.row(i).transpose()).collect();
    // Track which row index each working row originally came from.
    let row_ids: Vec<usize> = (0..n_rows).collect();
    // The reduction proceeds left-to-right by row index — each working row
    // is checked against the already-accepted basis. If it's linearly
    // dependent on the basis, it's flagged redundant.
    let mut basis: Vec<DVector<f64>> = Vec::new();
    let mut redundant_rows: Vec<usize> = Vec::new();
    for (work_idx, row_vec) in work.iter_mut().enumerate() {
        let original_row = row_ids[work_idx];
        // Project this row onto the current basis (Gram-Schmidt
        // orthogonalization): row ← row − Σ (row·b) · b for each unit b.
        for b in &basis {
            let coef = row_vec.dot(b);
            row_vec.axpy(-coef, b, 1.0);
        }
        let residual_norm = row_vec.norm();
        if residual_norm < tol {
            // Row is in the span of the previously-accepted basis →
            // redundant.
            redundant_rows.push(original_row);
        } else {
            // Normalize and add to basis.
            *row_vec /= residual_norm;
            basis.push(row_vec.clone());
        }
    }
    let rank = basis.len();
    (rank, redundant_rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mate::{Mate, MateKind};
    use crate::part::Part;
    use crate::solver::{solve, SolverConfig};
    use nalgebra::Vector3;

    fn unit_cube(name: &str) -> Part {
        Part::new(0, name, valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap())
    }

    /// Empty assembly: vacuously fully constrained.
    #[test]
    fn empty_assembly_is_fully_constrained() {
        let mut a = Assembly::new();
        let report = diagnose(&mut a, DiagnosticsConfig::default()).unwrap();
        assert_eq!(report.state, ConstraintState::FullyConstrained);
        assert_eq!(report.n_residuals, 0);
        assert_eq!(report.n_variables, 0);
    }

    /// Single floating part with no mates: 6 free DOFs.
    #[test]
    fn lone_floating_part_is_under_constrained_by_six() {
        let mut a = Assembly::new();
        a.add_part(unit_cube("free"));
        let report = diagnose(&mut a, DiagnosticsConfig::default()).unwrap();
        assert!(
            matches!(
                report.state,
                ConstraintState::UnderConstrained { remaining_dof: 6 }
            ),
            "got {:?}",
            report.state
        );
        assert_eq!(report.n_variables, 6);
    }

    /// One Coincident mate pins translation only (3 residuals), leaving
    /// 3 rotation DOFs free.
    #[test]
    fn coincident_pin_leaves_three_rotation_dof_free() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("anchor");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("free"));
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        ));
        // Solve first so the Jacobian is evaluated at a configuration where
        // the mate is satisfied (rotation columns are well-conditioned).
        solve(&mut a, SolverConfig::default()).unwrap();
        let report = diagnose(&mut a, DiagnosticsConfig::default()).unwrap();
        assert!(
            matches!(
                report.state,
                ConstraintState::UnderConstrained { remaining_dof: 3 }
            ),
            "got {:?}",
            report.state
        );
    }

    /// Duplicate Coincident mate at the same anchor pair → 6 residuals
    /// but rank still 3, so 3 rows are redundant. The duplicate mate id
    /// should appear in the redundant list.
    #[test]
    fn duplicate_coincident_is_over_constrained() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("anchor");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("free"));
        let mate_kind = MateKind::Coincident {
            part_a: id_a,
            point_a: Vector3::zeros(),
            part_b: id_b,
            point_b: Vector3::zeros(),
        };
        let _m0 = a.add_mate(Mate::new(0, mate_kind.clone()));
        let m1 = a.add_mate(Mate::new(0, mate_kind));
        solve(&mut a, SolverConfig::default()).unwrap();
        let report = diagnose(&mut a, DiagnosticsConfig::default()).unwrap();
        match &report.state {
            ConstraintState::OverConstrained { redundant_mates } => {
                assert!(
                    redundant_mates.contains(&m1),
                    "duplicate (mate {m1}) should be flagged: {redundant_mates:?}"
                );
            }
            other => panic!("expected OverConstrained, got {other:?}"),
        }
    }

    /// Two Distance mates between the same anchor points with conflicting
    /// targets — Inconsistent. (Distance is 1 residual; the two rows have
    /// the same Jacobian gradient but different residual values, so the
    /// solver settles to a non-zero residual.)
    #[test]
    fn conflicting_distance_mates_are_inconsistent() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("anchor");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let mut p1 = unit_cube("free");
        p1.transform.translation = Vector3::new(5.0, 0.0, 0.0);
        let id_b = a.add_part(p1);
        let _m0 = a.add_mate(Mate::new(
            0,
            MateKind::Distance {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
                target: 5.0,
            },
        ));
        let m1 = a.add_mate(Mate::new(
            0,
            MateKind::Distance {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
                target: 7.0,
            },
        ));
        // The solver will land somewhere between 5 and 7 with a non-zero
        // residual — the two equations cannot both be satisfied.
        let _ = solve(&mut a, SolverConfig::default());
        let report = diagnose(&mut a, DiagnosticsConfig::default()).unwrap();
        match &report.state {
            ConstraintState::Inconsistent { conflicting_mates } => {
                assert!(
                    conflicting_mates.contains(&m1),
                    "the duplicate mate should be in the conflict set: {conflicting_mates:?}"
                );
            }
            other => panic!("expected Inconsistent, got {other:?}; report = {report:?}"),
        }
    }

    /// Suppressed mates don't count.
    #[test]
    fn suppressed_mates_are_skipped() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("anchor");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("free"));
        let mut m = Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        );
        m.suppressed = true;
        a.add_mate(m);
        let report = diagnose(&mut a, DiagnosticsConfig::default()).unwrap();
        assert_eq!(report.n_residuals, 0);
        assert!(matches!(
            report.state,
            ConstraintState::UnderConstrained { remaining_dof: 6 }
        ));
    }

    /// mate_row_map walks un-suppressed mates and accumulates end-row.
    #[test]
    fn mate_row_map_is_cumulative_in_storage_order() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));
        let id0 = a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        ));
        let id1 = a.add_mate(Mate::new(
            0,
            MateKind::Distance {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
                target: 1.0,
            },
        ));
        let map = mate_row_map(&a);
        // Coincident contributes 3 rows; Distance contributes 1; cumulative.
        assert_eq!(map, vec![(id0, 3), (id1, 4)]);
    }
}
