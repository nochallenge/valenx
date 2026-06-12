//! Feature 16 — ensemble docking.
//!
//! A single crystal structure is one snapshot of a flexible protein.
//! *Ensemble docking* hedges against that: the same ligand is docked
//! against several receptor conformations (alternative crystal forms,
//! NMR models, MD snapshots) and the results are combined into one
//! per-ligand score.
//!
//! [`ensemble_dock`] runs the docking against each receptor
//! conformation and combines the per-conformation best scores by an
//! [`EnsembleMethod`]:
//!
//! - [`EnsembleMethod::Best`] — take the single most favourable score
//!   across all conformations (the conformation that binds best
//!   "wins"; the standard ensemble-docking rule).
//! - [`EnsembleMethod::Mean`] — average across conformations.
//! - [`EnsembleMethod::Boltzmann`] — a Boltzmann-weighted average,
//!   `Σ Eᵢ·wᵢ` with weights `wᵢ ∝ exp(−Eᵢ/kT)`, which down-weights
//!   high-energy conformations smoothly rather than discarding them.
//!   The result always lies between the single best score and the
//!   arithmetic mean.

use valenx_dock::ligand::Ligand;
use valenx_dock::receptor::Receptor;

use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;
use crate::search::driver::{rigid_dock, ScoredPose, SearchAlgorithm};

/// How per-conformation scores are combined into one ensemble score.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EnsembleMethod {
    /// The single best (most negative) score across all receptor
    /// conformations.
    Best,
    /// The arithmetic mean of the per-conformation best scores.
    Mean,
    /// A Boltzmann-weighted combination at temperature `kt`.
    Boltzmann {
        /// The `kT` value (energy units).
        kt: f64,
    },
}

impl EnsembleMethod {
    /// Combine a set of per-conformation best scores into one ensemble
    /// score. `scores` must be non-empty.
    fn combine(&self, scores: &[f64]) -> f64 {
        match *self {
            EnsembleMethod::Best => scores.iter().copied().fold(f64::INFINITY, f64::min),
            EnsembleMethod::Mean => scores.iter().sum::<f64>() / scores.len() as f64,
            EnsembleMethod::Boltzmann { kt } => {
                // Boltzmann-weighted average ⟨E⟩ = Σ Eᵢ·wᵢ with
                // wᵢ ∝ exp(−Eᵢ/kT). The exponent is shifted by the
                // minimum score for numerical stability (a common
                // constant factor cancels between numerator and the
                // partition sum). The weighted average always lies
                // between the minimum score and the arithmetic mean —
                // unlike the free energy −kT·ln Z, which can drop
                // below the minimum.
                let min = scores.iter().copied().fold(f64::INFINITY, f64::min);
                let mut z = 0.0;
                let mut weighted = 0.0;
                for &e in scores {
                    let w = (-(e - min) / kt).exp();
                    z += w;
                    weighted += e * w;
                }
                if z > 0.0 {
                    weighted / z
                } else {
                    min
                }
            }
        }
    }
}

/// The docking outcome against one receptor conformation.
#[derive(Clone, Debug)]
pub struct ConformationResult {
    /// Index of the receptor conformation this result belongs to.
    pub conformation: usize,
    /// The best pose found against this conformation.
    pub best_pose: ScoredPose,
}

/// The result of an ensemble-docking run.
#[derive(Clone, Debug)]
pub struct EnsembleResult {
    /// The combined ensemble score.
    pub ensemble_score: f64,
    /// Per-conformation docking results.
    pub per_conformation: Vec<ConformationResult>,
    /// Index of the conformation that produced the single best pose.
    pub best_conformation: usize,
}

impl EnsembleResult {
    /// The overall best pose across the whole ensemble.
    pub fn best_pose(&self) -> &ScoredPose {
        &self.per_conformation[self.best_conformation].best_pose
    }
}

/// Feature 16 — dock a ligand against an ensemble of receptor
/// conformations.
///
/// `receptors` is the conformation ensemble; each is docked
/// independently with [`rigid_dock`], and the per-conformation best
/// scores are combined by `method`.
///
/// Returns [`DockScreenError`] for an empty ensemble or an invalid
/// Boltzmann temperature.
pub fn ensemble_dock(
    receptors: &[Receptor],
    ligand: &Ligand,
    grid: &GridBox,
    method: EnsembleMethod,
    algorithm: SearchAlgorithm,
    runs_per_receptor: usize,
    seed: u64,
) -> Result<EnsembleResult> {
    if receptors.is_empty() {
        return Err(DockScreenError::invalid(
            "receptors",
            "ensemble must contain at least one receptor conformation",
        ));
    }
    if let EnsembleMethod::Boltzmann { kt } = method {
        if !kt.is_finite() || kt <= 0.0 {
            return Err(DockScreenError::invalid(
                "kt",
                "Boltzmann temperature must be positive",
            ));
        }
    }
    if runs_per_receptor == 0 {
        return Err(DockScreenError::invalid(
            "runs_per_receptor",
            "must dock against each conformation at least once",
        ));
    }

    let mut per_conformation: Vec<ConformationResult> = Vec::with_capacity(receptors.len());
    for (ci, receptor) in receptors.iter().enumerate() {
        let conf_seed = seed.wrapping_add((ci as u64).wrapping_mul(0x9E37_79B9));
        let run = rigid_dock(
            receptor,
            ligand,
            grid,
            algorithm,
            runs_per_receptor,
            conf_seed,
        )?;
        let best = run
            .poses
            .into_iter()
            .next()
            .ok_or_else(|| DockScreenError::invalid("search", "search returned no poses"))?;
        per_conformation.push(ConformationResult {
            conformation: ci,
            best_pose: best,
        });
    }

    let scores: Vec<f64> = per_conformation.iter().map(|c| c.best_pose.score).collect();
    let ensemble_score = method.combine(&scores);
    // The conformation with the single best (lowest) pose score.
    let best_conformation = per_conformation
        .iter()
        .enumerate()
        .min_by(|a, b| {
            a.1.best_pose
                .score
                .partial_cmp(&b.1.best_pose.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0);

    Ok(EnsembleResult {
        ensemble_score,
        per_conformation,
        best_conformation,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use valenx_dock::atom_type::Ad4AtomType;
    use valenx_dock::receptor::ReceptorAtom;

    fn carbon_receptor_at(x: f64) -> Receptor {
        Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::new(x, 0.0, 0.0),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        }
    }

    fn one_carbon_ligand() -> Ligand {
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
        Ligand::from_pdbqt(pdbqt).unwrap()
    }

    #[test]
    fn combine_best_takes_the_minimum() {
        assert_eq!(EnsembleMethod::Best.combine(&[-5.0, -8.0, -3.0]), -8.0);
    }

    #[test]
    fn combine_mean_averages() {
        assert!((EnsembleMethod::Mean.combine(&[-6.0, -8.0]) - -7.0).abs() < 1e-9);
    }

    #[test]
    fn combine_boltzmann_is_between_best_and_mean() {
        // The Boltzmann combination sits between the single best and
        // the arithmetic mean for any positive temperature.
        let scores = [-5.0, -8.0, -3.0];
        let best = EnsembleMethod::Best.combine(&scores);
        let mean = EnsembleMethod::Mean.combine(&scores);
        let boltz = EnsembleMethod::Boltzmann { kt: 1.0 }.combine(&scores);
        assert!(boltz <= mean + 1e-9 && boltz >= best - 1e-9);
    }

    #[test]
    fn rejects_empty_ensemble_and_bad_temperature() {
        let lig = one_carbon_ligand();
        let grid = GridBox::cubic([0.0; 3], 10.0).unwrap();
        assert!(ensemble_dock(
            &[],
            &lig,
            &grid,
            EnsembleMethod::Best,
            SearchAlgorithm::fast(),
            1,
            1
        )
        .is_err());
        let recs = vec![carbon_receptor_at(0.0)];
        assert!(ensemble_dock(
            &recs,
            &lig,
            &grid,
            EnsembleMethod::Boltzmann { kt: 0.0 },
            SearchAlgorithm::fast(),
            1,
            1
        )
        .is_err());
    }

    #[test]
    fn ensemble_docks_against_every_conformation() {
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([0.0; 3], [12.0; 3], 0.5).unwrap();
        let recs = vec![
            carbon_receptor_at(0.0),
            carbon_receptor_at(1.0),
            carbon_receptor_at(-1.0),
        ];
        let result = ensemble_dock(
            &recs,
            &lig,
            &grid,
            EnsembleMethod::Best,
            SearchAlgorithm::fast(),
            2,
            42,
        )
        .unwrap();
        assert_eq!(result.per_conformation.len(), 3);
        // ensemble_score (Best) equals the minimum per-conformation
        // score.
        let min = result
            .per_conformation
            .iter()
            .map(|c| c.best_pose.score)
            .fold(f64::INFINITY, f64::min);
        assert!((result.ensemble_score - min).abs() < 1e-9);
        assert!(result.best_pose().score.is_finite());
    }
}
