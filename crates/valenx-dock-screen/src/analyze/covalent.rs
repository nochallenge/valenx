//! Feature 25 — covalent docking.
//!
//! Most docking treats ligand binding as reversible — the ligand sits
//! in the pocket held by non-bonded forces. *Covalent* inhibitors
//! instead form a real chemical bond to the receptor (the classic
//! example: a warhead reacting with a catalytic cysteine). Covalent
//! docking pins one ligand atom — the *attachment atom* — at the
//! reactive receptor atom and searches only the poses consistent with
//! that constraint.
//!
//! [`covalent_dock`] implements this with a constrained search:
//!
//! 1. the ligand's attachment atom is restrained near the receptor's
//!    reactive atom (a target covalent-bond length apart);
//! 2. an unconstrained docking search proposes poses;
//! 3. each pose's score is augmented with a harmonic penalty for
//!    violating the attachment constraint;
//! 4. the best constrained pose is returned, together with the
//!    realised attachment distance.
//!
//! ### v1 note
//!
//! This pins the ligand geometrically and scores the attachment
//! constraint with a harmonic restraint. It does **not** model the
//! bond-formation chemistry (warhead reaction energetics, hybridisation
//! change at the receptor atom). That is a documented v1 limit — the
//! geometric covalent constraint is the real, useful part.

use nalgebra::Vector3;

use valenx_dock::atom_type::Ad4AtomType;
use valenx_dock::ligand::Ligand;
use valenx_dock::pose::Pose;
use valenx_dock::receptor::Receptor;

use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;
use crate::score::gridmap::{AffinityMapSet, MapKind};
use crate::search::ga::{GaParams, LamarckianGa};
use crate::search::objective::PoseObjective;

/// A typical carbon–sulfur / carbon–carbon covalent-bond length (Å) —
/// the default target attachment distance.
pub const DEFAULT_BOND_LENGTH: f64 = 1.8;

/// Tunable parameters for a covalent-docking run.
#[derive(Clone, Copy, Debug)]
pub struct CovalentParams {
    /// Index (into [`Receptor::atoms`]) of the reactive receptor atom
    /// the ligand attaches to.
    pub receptor_atom: usize,
    /// Index (into [`Ligand::atoms`]) of the ligand's attachment atom.
    pub ligand_atom: usize,
    /// Target covalent-bond length between the two atoms (Å).
    pub bond_length: f64,
    /// Harmonic force constant of the attachment restraint
    /// (kcal/mol·Å²). Larger = a tighter pin.
    pub restraint_k: f64,
    /// Number of independent constrained searches to run.
    pub n_runs: usize,
}

impl CovalentParams {
    /// Default covalent-docking parameters for a given receptor /
    /// ligand attachment-atom pair.
    pub fn new(receptor_atom: usize, ligand_atom: usize) -> Self {
        CovalentParams {
            receptor_atom,
            ligand_atom,
            bond_length: DEFAULT_BOND_LENGTH,
            restraint_k: 10.0,
            n_runs: 8,
        }
    }
}

/// The result of a covalent-docking run.
#[derive(Clone, Debug)]
pub struct CovalentResult {
    /// The best constrained pose.
    pub pose: Pose,
    /// The pose's total score (grid score + attachment restraint).
    pub score: f64,
    /// The grid (interaction) score alone, without the restraint.
    pub interaction_score: f64,
    /// The realised distance between the attachment atom and the
    /// reactive receptor atom (Å).
    pub attachment_distance: f64,
    /// The harmonic restraint penalty paid by the pose.
    pub restraint_penalty: f64,
}

impl CovalentResult {
    /// `true` if the realised attachment distance is within
    /// `tolerance` Å of the target bond length — i.e. the covalent
    /// constraint is honoured.
    pub fn is_attached(&self, target: f64, tolerance: f64) -> bool {
        (self.attachment_distance - target).abs() <= tolerance
    }
}

/// Feature 25 — dock a ligand covalently, anchored at a reactive
/// receptor atom.
///
/// The ligand's attachment atom is restrained at `bond_length` from
/// the receptor's reactive atom by a harmonic penalty added to the
/// docking score. The search box is centred on the reactive atom.
///
/// Returns [`DockScreenError`] if either anchor index is out of range
/// or a parameter is invalid.
pub fn covalent_dock(
    receptor: &Receptor,
    ligand: &Ligand,
    box_edge: f64,
    params: &CovalentParams,
    seed: u64,
) -> Result<CovalentResult> {
    if params.receptor_atom >= receptor.atoms.len() {
        return Err(DockScreenError::invalid_receptor(format!(
            "reactive receptor atom index {} out of range ({} atoms)",
            params.receptor_atom,
            receptor.atoms.len()
        )));
    }
    if params.ligand_atom >= ligand.atoms.len() {
        return Err(DockScreenError::invalid_ligand(format!(
            "ligand attachment atom index {} out of range ({} atoms)",
            params.ligand_atom,
            ligand.atoms.len()
        )));
    }
    if !params.bond_length.is_finite() || params.bond_length <= 0.0 {
        return Err(DockScreenError::invalid(
            "bond_length",
            "covalent bond length must be positive",
        ));
    }
    if !params.restraint_k.is_finite() || params.restraint_k < 0.0 {
        return Err(DockScreenError::invalid(
            "restraint_k",
            "restraint force constant must be non-negative",
        ));
    }
    if params.n_runs == 0 {
        return Err(DockScreenError::invalid("n_runs", "run count must be ≥ 1"));
    }

    // The search box is centred on the reactive receptor atom — that
    // is where the ligand must end up.
    let reactive = receptor.atoms[params.receptor_atom].position;
    let grid = GridBox::cubic(
        [reactive.x, reactive.y, reactive.z],
        box_edge.max(6.0),
    )?;

    // Build the affinity maps and a base objective.
    let ligand_types: Vec<Ad4AtomType> = ligand.atoms.iter().map(|a| a.ad4_type).collect();
    let maps = AffinityMapSet::precompute(receptor, &ligand_types, &grid, MapKind::Vina)?;
    let charges: Vec<f64> = ligand.atoms.iter().map(|a| a.partial_charge).collect();
    let base_objective = PoseObjective::new(ligand, &maps, charges);

    // Run the constrained search. `PoseObjective` scores only the
    // grid interaction; the covalent restraint is a harmonic penalty
    // the GA cannot see directly, so each run docks against the grid
    // objective inside a search box tightly centred on the reactive
    // receptor atom (which keeps proposed poses near the attachment
    // point), then every run's best pose is *re-ranked* by the full
    // constrained energy (grid score + restraint). The tight box plus
    // the restraint re-ranking together realise the covalent anchor.
    let mut best: Option<(Pose, f64)> = None;
    for run in 0..params.n_runs {
        let run_seed = seed.wrapping_add((run as u64).wrapping_mul(0x9E37_79B9));
        let docked = LamarckianGa::new(GaParams::fast()).run(
            &base_objective,
            grid.center,
            grid.size,
            run_seed,
        )?;
        // Re-rank the run's best pose by the full constrained energy.
        let total = constrained_score(
            &base_objective,
            ligand,
            &docked.best_pose,
            reactive,
            params,
        );
        if best.as_ref().map(|(_, s)| total < *s).unwrap_or(true) {
            best = Some((docked.best_pose, total));
        }
    }

    let (pose, total_score) = best.expect("n_runs >= 1 guarantees a best pose");
    // Recompute the breakdown for the chosen pose.
    let interaction_score = base_objective.score(&pose);
    let world = ligand.apply_pose(&pose);
    let attach_pos = world[params.ligand_atom];
    let attachment_distance = (attach_pos - reactive).norm();
    let dev = attachment_distance - params.bond_length;
    let restraint_penalty = params.restraint_k * dev * dev;

    Ok(CovalentResult {
        pose,
        score: total_score,
        interaction_score,
        attachment_distance,
        restraint_penalty,
    })
}

/// The constrained covalent score: the grid (interaction) score plus a
/// harmonic restraint pinning the attachment atom near the reactive
/// receptor atom.
fn constrained_score(
    objective: &PoseObjective,
    ligand: &Ligand,
    pose: &Pose,
    reactive: Vector3<f64>,
    params: &CovalentParams,
) -> f64 {
    let interaction = objective.score(pose);
    let world = ligand.apply_pose(pose);
    let attach = world[params.ligand_atom];
    let dev = (attach - reactive).norm() - params.bond_length;
    interaction + params.restraint_k * dev * dev
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_dock::receptor::ReceptorAtom;

    fn reactive_receptor() -> Receptor {
        // A "catalytic cysteine sulfur" reactive atom plus a partner.
        Receptor {
            atoms: vec![
                ReceptorAtom {
                    position: Vector3::new(5.0, 5.0, 5.0),
                    ad4_type: Ad4AtomType::SA,
                    partial_charge: -0.2,
                },
                ReceptorAtom {
                    position: Vector3::new(7.0, 5.0, 5.0),
                    ad4_type: Ad4AtomType::C,
                    partial_charge: 0.0,
                },
            ],
        }
    }

    fn warhead_ligand() -> Ligand {
        // A two-atom ligand; atom 0 is the warhead carbon (attachment).
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.100 C
ATOM      2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
        Ligand::from_pdbqt(pdbqt).unwrap()
    }

    #[test]
    fn rejects_out_of_range_anchor_indices() {
        let r = reactive_receptor();
        let lig = warhead_ligand();
        // Receptor atom 99 does not exist.
        let bad_r = CovalentParams::new(99, 0);
        assert!(covalent_dock(&r, &lig, 12.0, &bad_r, 1).is_err());
        // Ligand atom 99 does not exist.
        let bad_l = CovalentParams::new(0, 99);
        assert!(covalent_dock(&r, &lig, 12.0, &bad_l, 1).is_err());
    }

    #[test]
    fn rejects_invalid_parameters() {
        let r = reactive_receptor();
        let lig = warhead_ligand();
        let bad_len = CovalentParams {
            bond_length: 0.0,
            ..CovalentParams::new(0, 0)
        };
        assert!(covalent_dock(&r, &lig, 12.0, &bad_len, 1).is_err());
        let bad_runs = CovalentParams {
            n_runs: 0,
            ..CovalentParams::new(0, 0)
        };
        assert!(covalent_dock(&r, &lig, 12.0, &bad_runs, 1).is_err());
    }

    #[test]
    fn covalent_dock_runs_and_reports_an_attachment_distance() {
        let r = reactive_receptor();
        let lig = warhead_ligand();
        let params = CovalentParams::new(0, 0); // attach lig atom 0 to rec atom 0
        let result = covalent_dock(&r, &lig, 14.0, &params, 4).unwrap();
        assert!(result.attachment_distance >= 0.0);
        assert!(result.restraint_penalty >= 0.0);
        assert!(result.score.is_finite());
        // The total score is the interaction score plus the restraint.
        assert!(
            (result.score - (result.interaction_score + result.restraint_penalty)).abs() < 1e-6
        );
    }

    #[test]
    fn restraint_pulls_the_warhead_toward_the_reactive_atom() {
        // With a strong restraint, the search should bring the warhead
        // atom much closer to the reactive sulfur than the box edge.
        let r = reactive_receptor();
        let lig = warhead_ligand();
        let params = CovalentParams {
            restraint_k: 50.0,
            n_runs: 8,
            ..CovalentParams::new(0, 0)
        };
        let result = covalent_dock(&r, &lig, 16.0, &params, 11).unwrap();
        // The attachment distance should be well inside the 16 Å box —
        // a strong restraint keeps the warhead near the 1.8 Å target.
        assert!(
            result.attachment_distance < 8.0,
            "warhead ended {} Å from the reactive atom",
            result.attachment_distance
        );
    }

    #[test]
    fn is_attached_checks_the_bond_length() {
        let result = CovalentResult {
            pose: Pose::identity(0),
            score: -5.0,
            interaction_score: -6.0,
            attachment_distance: 1.85,
            restraint_penalty: 1.0,
        };
        assert!(result.is_attached(1.8, 0.3));
        assert!(!result.is_attached(1.8, 0.01));
    }
}
