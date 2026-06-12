//! Feature 21 — pose RMSD and redocking validation.
//!
//! How good *is* a docking protocol? The standard answer is a
//! *redocking* benchmark: take crystal complexes where the ligand's
//! true bound pose is known, dock each ligand back into its receptor,
//! and measure how often the top-ranked pose lands within 2 Å RMSD of
//! the crystal pose. That fraction — the **success rate** — is the
//! single number used to compare docking methods.
//!
//! This module provides:
//!
//! - [`pose_rmsd`] — heavy-atom RMSD between two poses of the same
//!   ligand (delegating to [`valenx_dock::cluster::rmsd`]);
//! - [`coordinate_rmsd`] — RMSD directly between two coordinate sets
//!   (for comparing a docked pose to crystal coordinates);
//! - [`redock_success_rate`] — run a redocking benchmark over a set of
//!   cases and report the success rate and per-case outcomes.
//!
//! The 2 Å threshold is the field convention; [`RedockBenchmark`]
//! lets the caller change it.

use nalgebra::Vector3;

use valenx_dock::cluster::rmsd as dock_rmsd;
use valenx_dock::ligand::Ligand;
use valenx_dock::pose::Pose;

use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;
use crate::search::driver::{rigid_dock, SearchAlgorithm};

/// The conventional redocking success-RMSD threshold (Å).
pub const SUCCESS_RMSD: f64 = 2.0;

/// Heavy-atom RMSD between two poses of the same ligand. Hydrogens are
/// excluded — the AutoDock convention. A thin wrapper over
/// [`valenx_dock::cluster::rmsd`].
pub fn pose_rmsd(ligand: &Ligand, a: &Pose, b: &Pose) -> f64 {
    dock_rmsd(ligand, a, b)
}

/// RMSD directly between two coordinate sets — element-wise, no
/// superposition (docking RMSD is computed in the receptor frame, so
/// the poses are already aligned).
///
/// Returns [`DockScreenError::Invalid`] if the two sets differ in
/// length or are empty.
pub fn coordinate_rmsd(a: &[Vector3<f64>], b: &[Vector3<f64>]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(DockScreenError::invalid(
            "coordinates",
            format!(
                "coordinate sets differ in length ({} vs {})",
                a.len(),
                b.len()
            ),
        ));
    }
    if a.is_empty() {
        return Err(DockScreenError::invalid(
            "coordinates",
            "cannot compute RMSD of empty coordinate sets",
        ));
    }
    let sum_sq: f64 = a
        .iter()
        .zip(b.iter())
        .map(|(p, q)| (p - q).norm_squared())
        .sum();
    Ok((sum_sq / a.len() as f64).sqrt())
}

/// One case of a redocking benchmark: a receptor, a ligand, the search
/// box, and the ligand's known reference (crystal) pose.
#[derive(Clone, Debug)]
pub struct RedockCase {
    /// A human-readable case identifier (e.g. a PDB id).
    pub name: String,
    /// Receptor PDBQT.
    pub receptor_pdbqt: String,
    /// Ligand PDBQT.
    pub ligand_pdbqt: String,
    /// The search box for this case.
    pub grid: GridBox,
    /// The known reference pose (the crystal binding mode).
    pub reference: Pose,
}

/// The redocking outcome for a single case.
#[derive(Clone, Debug)]
pub struct RedockOutcome {
    /// The case name.
    pub name: String,
    /// RMSD of the top-ranked docked pose to the reference pose (Å), or
    /// `None` if the case failed to dock.
    pub top_pose_rmsd: Option<f64>,
    /// `true` if `top_pose_rmsd` is below the success threshold.
    pub success: bool,
    /// A failure reason, or `None` if the case docked.
    pub failure: Option<String>,
}

/// The result of a redocking benchmark over several cases.
#[derive(Clone, Debug)]
pub struct RedockBenchmark {
    /// Per-case outcomes.
    pub outcomes: Vec<RedockOutcome>,
    /// Fraction `[0,1]` of cases whose top pose was within the
    /// threshold — the headline success rate.
    pub success_rate: f64,
    /// The RMSD threshold (Å) used.
    pub threshold: f64,
}

impl RedockBenchmark {
    /// Number of cases that docked successfully (regardless of RMSD).
    pub fn n_docked(&self) -> usize {
        self.outcomes.iter().filter(|o| o.failure.is_none()).count()
    }

    /// Mean top-pose RMSD over the cases that docked.
    pub fn mean_rmsd(&self) -> f64 {
        let rmsds: Vec<f64> = self
            .outcomes
            .iter()
            .filter_map(|o| o.top_pose_rmsd)
            .collect();
        if rmsds.is_empty() {
            return 0.0;
        }
        rmsds.iter().sum::<f64>() / rmsds.len() as f64
    }
}

/// Feature 21 — run a redocking benchmark.
///
/// Each case's ligand is docked back into its receptor; the top-ranked
/// pose is compared (by heavy-atom RMSD) to the case's reference pose.
/// A case "succeeds" when that RMSD is `< threshold`.
///
/// A case that fails to parse / dock is recorded as a non-success with
/// a `failure` reason rather than aborting the benchmark.
///
/// Returns [`DockScreenError::Invalid`] only for run-wide problems (an
/// empty case list, a non-positive threshold).
pub fn redock_success_rate(
    cases: &[RedockCase],
    threshold: f64,
    algorithm: SearchAlgorithm,
    runs_per_case: usize,
    seed: u64,
) -> Result<RedockBenchmark> {
    if cases.is_empty() {
        return Err(DockScreenError::invalid(
            "cases",
            "redocking benchmark needs at least one case",
        ));
    }
    if !threshold.is_finite() || threshold <= 0.0 {
        return Err(DockScreenError::invalid(
            "threshold",
            "success RMSD threshold must be positive",
        ));
    }
    if runs_per_case == 0 {
        return Err(DockScreenError::invalid(
            "runs_per_case",
            "must dock each case at least once",
        ));
    }

    let mut outcomes: Vec<RedockOutcome> = Vec::with_capacity(cases.len());
    for (i, case) in cases.iter().enumerate() {
        let case_seed = seed.wrapping_add((i as u64).wrapping_mul(0x9E37_79B9));
        outcomes.push(run_case(
            case,
            threshold,
            algorithm,
            runs_per_case,
            case_seed,
        ));
    }
    let n_success = outcomes.iter().filter(|o| o.success).count();
    let success_rate = n_success as f64 / cases.len() as f64;
    Ok(RedockBenchmark {
        outcomes,
        success_rate,
        threshold,
    })
}

/// Run a single redocking case, capturing any failure.
fn run_case(
    case: &RedockCase,
    threshold: f64,
    algorithm: SearchAlgorithm,
    runs: usize,
    seed: u64,
) -> RedockOutcome {
    let receptor = match valenx_dock::receptor::Receptor::from_pdbqt(&case.receptor_pdbqt) {
        Ok(r) => r,
        Err(e) => {
            return RedockOutcome {
                name: case.name.clone(),
                top_pose_rmsd: None,
                success: false,
                failure: Some(format!("receptor parse failed: {e}")),
            }
        }
    };
    let ligand = match Ligand::from_pdbqt(&case.ligand_pdbqt) {
        Ok(l) => l,
        Err(e) => {
            return RedockOutcome {
                name: case.name.clone(),
                top_pose_rmsd: None,
                success: false,
                failure: Some(format!("ligand parse failed: {e}")),
            }
        }
    };
    match rigid_dock(&receptor, &ligand, &case.grid, algorithm, runs, seed) {
        Ok(run) => match run.poses.first() {
            Some(top) => {
                let r = pose_rmsd(&ligand, &top.pose, &case.reference);
                RedockOutcome {
                    name: case.name.clone(),
                    top_pose_rmsd: Some(r),
                    success: r < threshold,
                    failure: None,
                }
            }
            None => RedockOutcome {
                name: case.name.clone(),
                top_pose_rmsd: None,
                success: false,
                failure: Some("search returned no poses".into()),
            },
        },
        Err(e) => RedockOutcome {
            name: case.name.clone(),
            top_pose_rmsd: None,
            success: false,
            failure: Some(format!("dock failed: {e}")),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECEPTOR: &str =
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  CB  ALA A   1       2.000   0.000   0.000  1.00  0.00     0.000 C
";
    const LIGAND: &str = "ROOT
ATOM      1  C1  LIG A   1       1.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";

    fn one_atom_ligand() -> Ligand {
        Ligand::from_pdbqt(LIGAND).unwrap()
    }

    #[test]
    fn pose_rmsd_of_identical_poses_is_zero() {
        let lig = one_atom_ligand();
        let p = Pose::identity(0);
        assert_eq!(pose_rmsd(&lig, &p, &p), 0.0);
    }

    #[test]
    fn pose_rmsd_pure_translation_equals_distance() {
        let lig = one_atom_ligand();
        let mut a = Pose::identity(0);
        let mut b = Pose::identity(0);
        a.translation = Vector3::zeros();
        b.translation = Vector3::new(3.0, 4.0, 0.0);
        assert!((pose_rmsd(&lig, &a, &b) - 5.0).abs() < 1e-9);
    }

    #[test]
    fn coordinate_rmsd_basic() {
        let a = vec![Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)];
        let b = vec![Vector3::new(3.0, 0.0, 0.0), Vector3::new(4.0, 0.0, 0.0)];
        // Each pair is 3 Å apart → RMSD 3.
        assert!((coordinate_rmsd(&a, &b).unwrap() - 3.0).abs() < 1e-9);
    }

    #[test]
    fn coordinate_rmsd_rejects_mismatch_and_empty() {
        assert!(coordinate_rmsd(&[Vector3::zeros()], &[]).is_err());
        assert!(coordinate_rmsd(&[], &[]).is_err());
    }

    #[test]
    fn redock_rejects_degenerate_inputs() {
        assert!(redock_success_rate(&[], 2.0, SearchAlgorithm::fast(), 1, 1).is_err());
        let case = RedockCase {
            name: "x".into(),
            receptor_pdbqt: RECEPTOR.into(),
            ligand_pdbqt: LIGAND.into(),
            grid: GridBox::cubic([1.0, 0.0, 0.0], 8.0).unwrap(),
            reference: Pose::identity(0),
        };
        // Bad threshold.
        assert!(redock_success_rate(
            std::slice::from_ref(&case),
            0.0,
            SearchAlgorithm::fast(),
            1,
            1
        )
        .is_err());
    }

    #[test]
    fn redock_benchmark_runs_and_reports_a_rate() {
        let case = RedockCase {
            name: "test-complex".into(),
            receptor_pdbqt: RECEPTOR.into(),
            ligand_pdbqt: LIGAND.into(),
            grid: GridBox::with_spacing([1.0, 0.0, 0.0], [12.0; 3], 0.5).unwrap(),
            reference: {
                let mut p = Pose::identity(0);
                p.translation = Vector3::new(0.5, 0.0, 0.0);
                p
            },
        };
        let bench =
            redock_success_rate(&[case], SUCCESS_RMSD, SearchAlgorithm::fast(), 3, 42).unwrap();
        assert_eq!(bench.outcomes.len(), 1);
        assert!((0.0..=1.0).contains(&bench.success_rate));
        assert_eq!(bench.threshold, SUCCESS_RMSD);
        assert_eq!(bench.n_docked(), 1);
        assert!(bench.mean_rmsd() >= 0.0);
    }

    #[test]
    fn a_failed_case_is_recorded_not_propagated() {
        let bad_case = RedockCase {
            name: "broken".into(),
            receptor_pdbqt: "not a receptor".into(),
            ligand_pdbqt: LIGAND.into(),
            grid: GridBox::cubic([0.0; 3], 8.0).unwrap(),
            reference: Pose::identity(0),
        };
        let bench = redock_success_rate(&[bad_case], 2.0, SearchAlgorithm::fast(), 1, 1).unwrap();
        assert_eq!(bench.outcomes.len(), 1);
        assert!(!bench.outcomes[0].success);
        assert!(bench.outcomes[0].failure.is_some());
        assert_eq!(bench.n_docked(), 0);
    }
}
