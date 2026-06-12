//! Feature 14 — batch virtual screening over a ligand library.
//!
//! Virtual screening docks an entire compound library against one
//! receptor and ranks the hits. The mechanics are simple — dock each
//! ligand, keep its best score — but two things matter for a real
//! screen:
//!
//! - **throughput** — ligands are independent, so the library is
//!   docked in parallel across CPU cores with `rayon`;
//! - **robustness** — a single malformed ligand must not abort the
//!   whole run, so a per-ligand parse / dock failure is captured as a
//!   failed [`ScreenEntry`] rather than propagated.
//!
//! The output is a [`ScreenResult`]: every ligand's best score, sorted
//! best-first, plus the count of ligands that failed.

use rayon::prelude::*;

use valenx_dock::ligand::Ligand;
use valenx_dock::receptor::Receptor;

use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;
use crate::search::driver::{rigid_dock, ScoredPose, SearchAlgorithm};

/// One ligand in a screening library: a name and its PDBQT text.
#[derive(Clone, Debug)]
pub struct LibraryEntry {
    /// A human-readable identifier (compound id, catalogue number, …).
    pub name: String,
    /// The ligand's PDBQT document.
    pub pdbqt: String,
}

impl LibraryEntry {
    /// Build a library entry.
    pub fn new(name: impl Into<String>, pdbqt: impl Into<String>) -> Self {
        LibraryEntry {
            name: name.into(),
            pdbqt: pdbqt.into(),
        }
    }
}

/// The screening outcome for a single ligand.
#[derive(Clone, Debug)]
pub struct ScreenEntry {
    /// The ligand's name (from the [`LibraryEntry`]).
    pub name: String,
    /// The best docking score, or `None` if the ligand failed to
    /// parse / dock.
    pub best_score: Option<f64>,
    /// The best pose, or `None` on failure.
    pub best_pose: Option<ScoredPose>,
    /// Number of distinct poses the search returned (`0` on failure).
    pub n_poses: usize,
    /// A failure reason, or `None` if the ligand docked successfully.
    pub failure: Option<String>,
}

impl ScreenEntry {
    /// `true` if the ligand docked without error.
    pub fn succeeded(&self) -> bool {
        self.failure.is_none()
    }
}

/// The result of screening a whole library.
#[derive(Clone, Debug)]
pub struct ScreenResult {
    /// Per-ligand outcomes, sorted by best score (best first; failed
    /// ligands sort last).
    pub entries: Vec<ScreenEntry>,
    /// Number of ligands that failed to parse / dock.
    pub n_failed: usize,
}

impl ScreenResult {
    /// The top `n` successful entries by score.
    pub fn top(&self, n: usize) -> Vec<&ScreenEntry> {
        self.entries
            .iter()
            .filter(|e| e.succeeded())
            .take(n)
            .collect()
    }

    /// Number of ligands screened.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// `true` if no ligands were screened.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// The single best-scoring ligand, if any succeeded.
    pub fn best(&self) -> Option<&ScreenEntry> {
        self.entries.iter().find(|e| e.succeeded())
    }
}

/// Feature 14 — screen a ligand library against a receptor.
///
/// Every ligand is docked with [`rigid_dock`]; the library is
/// processed in parallel. A ligand that fails to parse or dock yields
/// a failed [`ScreenEntry`] (with a `failure` reason) rather than
/// aborting the screen. The result is sorted best-score-first.
///
/// `runs_per_ligand` is how many independent searches each ligand
/// gets (more runs → better global-minimum coverage, slower).
///
/// Returns [`DockScreenError`] only for run-wide problems (an empty
/// library, an empty receptor); per-ligand failures are reported
/// inside the [`ScreenResult`].
pub fn screen_library(
    receptor: &Receptor,
    library: &[LibraryEntry],
    grid: &GridBox,
    algorithm: SearchAlgorithm,
    runs_per_ligand: usize,
    seed: u64,
) -> Result<ScreenResult> {
    if receptor.atoms.is_empty() {
        return Err(DockScreenError::invalid_receptor("receptor has no atoms"));
    }
    if library.is_empty() {
        return Err(DockScreenError::invalid(
            "library",
            "cannot screen an empty ligand library",
        ));
    }
    if runs_per_ligand == 0 {
        return Err(DockScreenError::invalid(
            "runs_per_ligand",
            "must dock each ligand at least once",
        ));
    }

    let mut entries: Vec<ScreenEntry> = library
        .par_iter()
        .enumerate()
        .map(|(i, entry)| {
            let ligand_seed = seed.wrapping_add((i as u64).wrapping_mul(0x9E37_79B9));
            dock_one(
                receptor,
                entry,
                grid,
                algorithm,
                runs_per_ligand,
                ligand_seed,
            )
        })
        .collect();

    let n_failed = entries.iter().filter(|e| !e.succeeded()).count();
    // Sort: successes by ascending score, failures last.
    entries.sort_by(|a, b| match (a.best_score, b.best_score) {
        (Some(sa), Some(sb)) => sa.partial_cmp(&sb).unwrap_or(std::cmp::Ordering::Equal),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });
    Ok(ScreenResult { entries, n_failed })
}

/// Dock one library entry, capturing any failure as a failed
/// [`ScreenEntry`] instead of propagating it.
fn dock_one(
    receptor: &Receptor,
    entry: &LibraryEntry,
    grid: &GridBox,
    algorithm: SearchAlgorithm,
    runs: usize,
    seed: u64,
) -> ScreenEntry {
    let ligand = match Ligand::from_pdbqt(&entry.pdbqt) {
        Ok(l) => l,
        Err(e) => {
            return ScreenEntry {
                name: entry.name.clone(),
                best_score: None,
                best_pose: None,
                n_poses: 0,
                failure: Some(format!("ligand parse failed: {e}")),
            }
        }
    };
    match rigid_dock(receptor, &ligand, grid, algorithm, runs, seed) {
        Ok(run) => {
            let best = run.poses.first().cloned();
            ScreenEntry {
                name: entry.name.clone(),
                best_score: best.as_ref().map(|p| p.score),
                best_pose: best,
                n_poses: run.poses.len(),
                failure: None,
            }
        }
        Err(e) => ScreenEntry {
            name: entry.name.clone(),
            best_score: None,
            best_pose: None,
            n_poses: 0,
            failure: Some(format!("dock failed: {e}")),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use valenx_dock::atom_type::Ad4AtomType;
    use valenx_dock::receptor::ReceptorAtom;

    fn carbon_receptor() -> Receptor {
        Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        }
    }

    const GOOD_LIGAND: &str = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";

    #[test]
    fn rejects_empty_library_and_receptor() {
        let r = carbon_receptor();
        let grid = GridBox::cubic([0.0; 3], 10.0).unwrap();
        assert!(screen_library(&r, &[], &grid, SearchAlgorithm::fast(), 1, 1).is_err());
        let lib = vec![LibraryEntry::new("x", GOOD_LIGAND)];
        assert!(screen_library(
            &Receptor::default(),
            &lib,
            &grid,
            SearchAlgorithm::fast(),
            1,
            1
        )
        .is_err());
    }

    #[test]
    fn screens_a_small_library_and_ranks_it() {
        let r = carbon_receptor();
        let grid = GridBox::with_spacing([0.0; 3], [12.0; 3], 0.5).unwrap();
        let lib = vec![
            LibraryEntry::new("cpd-1", GOOD_LIGAND),
            LibraryEntry::new("cpd-2", GOOD_LIGAND),
            LibraryEntry::new("cpd-3", GOOD_LIGAND),
        ];
        let result = screen_library(&r, &lib, &grid, SearchAlgorithm::fast(), 2, 42).unwrap();
        assert_eq!(result.len(), 3);
        assert_eq!(result.n_failed, 0);
        // Entries sorted ascending by score.
        let scores: Vec<f64> = result.entries.iter().filter_map(|e| e.best_score).collect();
        for w in scores.windows(2) {
            assert!(w[0] <= w[1], "screen result not sorted: {scores:?}");
        }
        assert!(result.best().is_some());
    }

    #[test]
    fn a_malformed_ligand_does_not_abort_the_screen() {
        let r = carbon_receptor();
        let grid = GridBox::with_spacing([0.0; 3], [12.0; 3], 0.5).unwrap();
        let lib = vec![
            LibraryEntry::new("good", GOOD_LIGAND),
            // No ROOT record → parse failure.
            LibraryEntry::new("bad", "TORSDOF 0\n"),
        ];
        let result = screen_library(&r, &lib, &grid, SearchAlgorithm::fast(), 1, 1).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result.n_failed, 1);
        // The good one succeeded and sorts first.
        assert!(result.entries[0].succeeded());
        // The failed one carries a reason and sorts last.
        let failed = result.entries.iter().find(|e| !e.succeeded()).unwrap();
        assert!(failed.failure.as_ref().unwrap().contains("parse"));
    }

    #[test]
    fn top_returns_only_successes() {
        let r = carbon_receptor();
        let grid = GridBox::with_spacing([0.0; 3], [12.0; 3], 0.5).unwrap();
        let lib = vec![
            LibraryEntry::new("good", GOOD_LIGAND),
            LibraryEntry::new("bad", "garbage\n"),
        ];
        let result = screen_library(&r, &lib, &grid, SearchAlgorithm::fast(), 1, 1).unwrap();
        let top = result.top(10);
        assert!(top.iter().all(|e| e.succeeded()));
    }
}
