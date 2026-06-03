//! Feature 30 — top-level docking and screening drivers.
//!
//! This is the crate's front door: two one-call entry points that take
//! PDBQT text in and return a structured report out, plus a registry
//! that probes which external neural-network tools are installed.
//!
//! - [`dock`] — dock a single ligand against a receptor and return a
//!   [`DockingReport`] (clustered, ranked poses + an interaction
//!   fingerprint for the top pose).
//! - [`screen`] — dock a whole ligand library and return a
//!   [`ScreeningReport`] (per-ligand best scores, ranked).
//! - [`AdapterRegistry`] — probe the host `PATH` and report which of
//!   the structure-prediction / generative-design / NN-docking /
//!   cryo-EM external tools are available.
//!
//! The drivers wire together the preparation, scoring, search,
//! clustering and analysis layers; [`DockParams`] bundles the few
//! knobs a caller usually wants to set.

use valenx_dock::ligand::Ligand;
use valenx_dock::receptor::Receptor;

use crate::adapters::common::{find_executable, ToolStatus};
use crate::adapters::cryo_em::CryoEmTool;
use crate::adapters::generative_design::GenerativeTool;
use crate::adapters::nn_docking::NnDockingTool;
use crate::adapters::structure_prediction::StructurePredictionTool;
use crate::analyze::fingerprint::{interaction_fingerprint, InteractionFingerprint};
use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;
use crate::screen::batch::{screen_library, LibraryEntry, ScreenEntry};
use crate::screen::cluster::{cluster_poses, PoseCluster};
use crate::search::driver::{rigid_dock, ScoredPose, SearchAlgorithm};

/// Tunable parameters for the top-level [`dock`] / [`screen`] drivers.
#[derive(Clone, Copy, Debug)]
pub struct DockParams {
    /// The search box. If `None`, a box is derived to enclose the
    /// receptor (blind docking).
    pub grid: Option<GridBox>,
    /// Which search algorithm to run.
    pub algorithm: SearchAlgorithm,
    /// Number of independent searches per ligand.
    pub n_runs: usize,
    /// RMSD cutoff (Å) for clustering the returned poses.
    pub cluster_rmsd: f64,
    /// Reproducibility seed.
    pub seed: u64,
}

impl Default for DockParams {
    fn default() -> Self {
        DockParams {
            grid: None,
            algorithm: SearchAlgorithm::fast(),
            n_runs: 8,
            cluster_rmsd: 2.0,
            seed: 0,
        }
    }
}

impl DockParams {
    /// A small, fast configuration for previews and unit tests.
    pub fn fast() -> Self {
        DockParams {
            grid: None,
            algorithm: SearchAlgorithm::fast(),
            n_runs: 3,
            cluster_rmsd: 2.0,
            seed: 0,
        }
    }

    /// A configuration with an explicit search box.
    pub fn with_grid(mut self, grid: GridBox) -> Self {
        self.grid = Some(grid);
        self
    }
}

/// The result of a single top-level docking run.
#[derive(Clone, Debug)]
pub struct DockingReport {
    /// Distinct binding modes (clusters), ranked best-first.
    pub clusters: Vec<PoseCluster>,
    /// The flat list of all ranked poses (best-first).
    pub poses: Vec<ScoredPose>,
    /// The interaction fingerprint of the top-ranked pose.
    pub top_fingerprint: InteractionFingerprint,
    /// The search box that was used.
    pub grid: GridBox,
}

impl DockingReport {
    /// The single best pose, if any.
    pub fn best(&self) -> Option<&ScoredPose> {
        self.poses.first()
    }

    /// The best docking score, if any.
    pub fn best_score(&self) -> Option<f64> {
        self.poses.first().map(|p| p.score)
    }

    /// Number of distinct binding modes found.
    pub fn n_binding_modes(&self) -> usize {
        self.clusters.len()
    }
}

/// Feature 30 — dock a single ligand against a receptor.
///
/// Parses the receptor and ligand PDBQT, picks (or derives) a search
/// box, runs the search, clusters the poses, and fingerprints the top
/// pose.
///
/// Returns [`DockScreenError`] if the inputs do not parse or a
/// parameter is invalid.
pub fn dock(
    receptor_pdbqt: &str,
    ligand_pdbqt: &str,
    params: &DockParams,
) -> Result<DockingReport> {
    let receptor = Receptor::from_pdbqt(receptor_pdbqt)?;
    let ligand = Ligand::from_pdbqt(ligand_pdbqt)?;

    // Pick the search box: the caller's, or one enclosing the receptor.
    let grid = match params.grid {
        Some(g) => g,
        None => {
            let pts: Vec<nalgebra::Vector3<f64>> =
                receptor.atoms.iter().map(|a| a.position).collect();
            GridBox::enclosing(&pts, 4.0)?
        }
    };

    let run = rigid_dock(
        &receptor,
        &ligand,
        &grid,
        params.algorithm,
        params.n_runs,
        params.seed,
    )?;
    let clusters = cluster_poses(&ligand, &run.poses, params.cluster_rmsd)?;

    // Fingerprint the top pose.
    let top_fingerprint = match run.poses.first() {
        Some(top) => {
            let world = ligand.apply_pose(&top.pose);
            let atoms: Vec<_> = world
                .iter()
                .zip(ligand.atoms.iter())
                .map(|(p, a)| (*p, a.ad4_type, a.partial_charge))
                .collect();
            interaction_fingerprint(&receptor, &atoms)
        }
        None => InteractionFingerprint::default(),
    };

    Ok(DockingReport {
        clusters,
        poses: run.poses,
        top_fingerprint,
        grid,
    })
}

/// One ligand's entry in a screening report.
#[derive(Clone, Debug)]
pub struct ScreeningHit {
    /// The ligand's name.
    pub name: String,
    /// The best docking score, or `None` if the ligand failed.
    pub best_score: Option<f64>,
    /// A failure reason, or `None` if the ligand docked.
    pub failure: Option<String>,
}

impl From<&ScreenEntry> for ScreeningHit {
    fn from(e: &ScreenEntry) -> Self {
        ScreeningHit {
            name: e.name.clone(),
            best_score: e.best_score,
            failure: e.failure.clone(),
        }
    }
}

/// The result of a top-level virtual-screening run.
#[derive(Clone, Debug)]
pub struct ScreeningReport {
    /// Per-ligand hits, ranked best-score-first (failures last).
    pub hits: Vec<ScreeningHit>,
    /// Number of ligands screened.
    pub n_screened: usize,
    /// Number of ligands that failed to parse / dock.
    pub n_failed: usize,
}

impl ScreeningReport {
    /// The top `n` successful hits.
    pub fn top(&self, n: usize) -> Vec<&ScreeningHit> {
        self.hits
            .iter()
            .filter(|h| h.failure.is_none())
            .take(n)
            .collect()
    }

    /// The single best-scoring ligand, if any succeeded.
    pub fn best(&self) -> Option<&ScreeningHit> {
        self.hits.iter().find(|h| h.failure.is_none())
    }
}

/// Feature 30 — screen a ligand library against a receptor.
///
/// `library` is `(name, ligand_pdbqt)` pairs. Each ligand is docked
/// and the library is ranked best-score-first. A ligand that fails to
/// parse / dock is recorded as a failed hit rather than aborting.
///
/// Returns [`DockScreenError`] only for run-wide problems (a bad
/// receptor, an empty library).
pub fn screen(
    receptor_pdbqt: &str,
    library: &[(String, String)],
    params: &DockParams,
) -> Result<ScreeningReport> {
    let receptor = Receptor::from_pdbqt(receptor_pdbqt)?;
    if library.is_empty() {
        return Err(DockScreenError::invalid(
            "library",
            "cannot screen an empty ligand library",
        ));
    }
    let grid = match params.grid {
        Some(g) => g,
        None => {
            let pts: Vec<nalgebra::Vector3<f64>> =
                receptor.atoms.iter().map(|a| a.position).collect();
            GridBox::enclosing(&pts, 4.0)?
        }
    };
    let entries: Vec<LibraryEntry> = library
        .iter()
        .map(|(name, pdbqt)| LibraryEntry::new(name.clone(), pdbqt.clone()))
        .collect();
    let result = screen_library(
        &receptor,
        &entries,
        &grid,
        params.algorithm,
        params.n_runs,
        params.seed,
    )?;
    let hits: Vec<ScreeningHit> = result.entries.iter().map(ScreeningHit::from).collect();
    Ok(ScreeningReport {
        n_screened: hits.len(),
        n_failed: result.n_failed,
        hits,
    })
}

/// The availability of one external neural-network tool.
#[derive(Clone, Debug, PartialEq)]
pub struct ToolAvailability {
    /// A stable tool label.
    pub label: &'static str,
    /// A human-readable display name.
    pub display_name: &'static str,
    /// Whether the tool's binary was found on `PATH`.
    pub status: ToolStatus,
}

impl ToolAvailability {
    /// `true` if the tool is installed.
    pub fn is_available(&self) -> bool {
        self.status.is_available()
    }
}

/// Feature 30 — a registry that probes which external neural-network
/// tools are installed on the host.
///
/// The classical docking in this crate needs no external tools — it is
/// a real native implementation. The *neural-network* tools
/// ([`crate::adapters`]) do, and a caller (or the Valenx UI) wants to
/// know up-front which are present. [`AdapterRegistry::probe`] checks
/// `PATH` for every adapter-wrapped tool and returns one
/// [`ToolAvailability`] each.
#[derive(Clone, Debug, Default)]
pub struct AdapterRegistry {
    /// Structure-prediction tool availabilities.
    pub structure_prediction: Vec<ToolAvailability>,
    /// Generative-design tool availabilities.
    pub generative_design: Vec<ToolAvailability>,
    /// Neural-network docking tool availabilities.
    pub nn_docking: Vec<ToolAvailability>,
    /// Cryo-EM reconstruction tool availabilities.
    pub cryo_em: Vec<ToolAvailability>,
}

impl AdapterRegistry {
    /// Probe the host `PATH` for every adapter-wrapped external tool.
    ///
    /// This does not run anything — it only checks whether each tool's
    /// binary exists. Safe to call cheaply at startup.
    pub fn probe() -> Self {
        let sp = [
            StructurePredictionTool::AlphaFold2,
            StructurePredictionTool::AlphaFold3,
            StructurePredictionTool::ColabFold,
            StructurePredictionTool::EsmFold,
            StructurePredictionTool::RoseTTAFold,
            StructurePredictionTool::OmegaFold,
            StructurePredictionTool::Boltz,
            StructurePredictionTool::Chai,
        ];
        let gd = [
            GenerativeTool::RfDiffusion,
            GenerativeTool::ProteinMpnn,
            GenerativeTool::EsmIf,
            GenerativeTool::Chroma,
        ];
        let nd = [NnDockingTool::DiffDock, NnDockingTool::Gnina];
        let ce = [CryoEmTool::Relion, CryoEmTool::CryoSparc, CryoEmTool::Eman2];

        AdapterRegistry {
            structure_prediction: sp
                .iter()
                .map(|t| probe_structure_prediction(*t))
                .collect(),
            generative_design: gd.iter().map(|t| probe_generative(*t)).collect(),
            nn_docking: nd.iter().map(|t| probe_nn_docking(*t)).collect(),
            cryo_em: ce.iter().map(|t| probe_cryo_em(*t)).collect(),
        }
    }

    /// Every probed tool, across all four categories.
    pub fn all(&self) -> Vec<&ToolAvailability> {
        self.structure_prediction
            .iter()
            .chain(self.generative_design.iter())
            .chain(self.nn_docking.iter())
            .chain(self.cryo_em.iter())
            .collect()
    }

    /// The total number of external tools probed.
    pub fn total_tools(&self) -> usize {
        self.all().len()
    }

    /// The number of external tools found installed on the host.
    pub fn available_count(&self) -> usize {
        self.all().iter().filter(|t| t.is_available()).count()
    }
}

// --- per-category probe helpers --------------------------------------
// Each re-uses the same `find_executable` PATH probe the adapters use.

fn probe_structure_prediction(tool: StructurePredictionTool) -> ToolAvailability {
    let candidates: &[&str] = match tool {
        StructurePredictionTool::AlphaFold2 => &["run_alphafold.py", "alphafold"],
        StructurePredictionTool::AlphaFold3 => &["run_alphafold", "alphafold3"],
        StructurePredictionTool::ColabFold => &["colabfold_batch"],
        StructurePredictionTool::EsmFold => &["esm-fold", "esmfold"],
        StructurePredictionTool::RoseTTAFold => &["rosettafold", "run_RF2.sh"],
        StructurePredictionTool::OmegaFold => &["omegafold"],
        StructurePredictionTool::Boltz => &["boltz"],
        StructurePredictionTool::Chai => &["chai", "chai-lab"],
    };
    ToolAvailability {
        label: tool.label(),
        display_name: tool.display_name(),
        status: find_executable(candidates),
    }
}

fn probe_generative(tool: GenerativeTool) -> ToolAvailability {
    let candidates: &[&str] = match tool {
        GenerativeTool::RfDiffusion => &["run_inference.py", "rfdiffusion"],
        GenerativeTool::ProteinMpnn => &["protein_mpnn_run.py", "proteinmpnn"],
        GenerativeTool::EsmIf => &["esm-if", "esm_inverse_folding"],
        GenerativeTool::Chroma => &["chroma"],
    };
    ToolAvailability {
        label: tool.label(),
        display_name: tool.display_name(),
        status: find_executable(candidates),
    }
}

fn probe_nn_docking(tool: NnDockingTool) -> ToolAvailability {
    let candidates: &[&str] = match tool {
        NnDockingTool::DiffDock => &["diffdock", "run_diffdock.py"],
        NnDockingTool::Gnina => &["gnina"],
    };
    ToolAvailability {
        label: tool.label(),
        display_name: tool.display_name(),
        status: find_executable(candidates),
    }
}

fn probe_cryo_em(tool: CryoEmTool) -> ToolAvailability {
    let candidates: &[&str] = match tool {
        CryoEmTool::Relion => &["relion_refine", "relion_refine_mpi"],
        CryoEmTool::CryoSparc => &["cryosparcm"],
        CryoEmTool::Eman2 => &["e2refine_easy.py", "e2refine2d.py"],
    };
    ToolAvailability {
        label: tool.label(),
        display_name: tool.display_name(),
        status: find_executable(candidates),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECEPTOR: &str = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  CB  ALA A   1       2.000   0.000   0.000  1.00  0.00     0.000 C
";
    const LIGAND: &str = "ROOT
ATOM      1  C1  LIG A   1       1.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";

    #[test]
    fn dock_returns_a_report_with_poses_and_a_fingerprint() {
        let report = dock(RECEPTOR, LIGAND, &DockParams::fast()).unwrap();
        assert!(!report.poses.is_empty());
        assert!(report.best().is_some());
        assert!(report.best_score().is_some());
        assert!(report.n_binding_modes() >= 1);
        // The fingerprint exists (may be empty for this tiny case).
        let _ = report.top_fingerprint.total();
    }

    #[test]
    fn dock_with_an_explicit_grid_uses_it() {
        let grid = GridBox::with_spacing([1.0, 0.0, 0.0], [12.0; 3], 0.5).unwrap();
        let params = DockParams::fast().with_grid(grid);
        let report = dock(RECEPTOR, LIGAND, &params).unwrap();
        assert_eq!(report.grid.center, grid.center);
    }

    #[test]
    fn dock_rejects_bad_input() {
        // A receptor that does not parse.
        assert!(dock("not a receptor", LIGAND, &DockParams::fast()).is_err());
        // A ligand with no ROOT.
        assert!(dock(RECEPTOR, "TORSDOF 0\n", &DockParams::fast()).is_err());
    }

    #[test]
    fn screen_ranks_a_small_library() {
        let library = vec![
            ("cpd-a".to_string(), LIGAND.to_string()),
            ("cpd-b".to_string(), LIGAND.to_string()),
        ];
        let report = screen(RECEPTOR, &library, &DockParams::fast()).unwrap();
        assert_eq!(report.n_screened, 2);
        assert_eq!(report.n_failed, 0);
        assert!(report.best().is_some());
        // Successful hits are sorted ascending by score.
        let scores: Vec<f64> = report.hits.iter().filter_map(|h| h.best_score).collect();
        for w in scores.windows(2) {
            assert!(w[0] <= w[1]);
        }
    }

    #[test]
    fn screen_rejects_an_empty_library() {
        assert!(screen(RECEPTOR, &[], &DockParams::fast()).is_err());
    }

    #[test]
    fn screen_records_a_failed_ligand() {
        let library = vec![
            ("good".to_string(), LIGAND.to_string()),
            ("bad".to_string(), "garbage\n".to_string()),
        ];
        let report = screen(RECEPTOR, &library, &DockParams::fast()).unwrap();
        assert_eq!(report.n_failed, 1);
        assert!(report.top(10).iter().all(|h| h.failure.is_none()));
    }

    #[test]
    fn adapter_registry_probes_every_tool() {
        let reg = AdapterRegistry::probe();
        // 8 structure-prediction + 4 generative + 2 NN-docking + 3
        // cryo-EM = 17 tools probed.
        assert_eq!(reg.structure_prediction.len(), 8);
        assert_eq!(reg.generative_design.len(), 4);
        assert_eq!(reg.nn_docking.len(), 2);
        assert_eq!(reg.cryo_em.len(), 3);
        assert_eq!(reg.total_tools(), 17);
        // available_count is between 0 and total — on a dev / CI host
        // it is almost always 0, which is correct and honest.
        assert!(reg.available_count() <= reg.total_tools());
    }

    #[test]
    fn adapter_registry_all_concatenates_categories() {
        let reg = AdapterRegistry::probe();
        assert_eq!(reg.all().len(), reg.total_tools());
        // Every probed tool carries a non-empty label and name.
        for t in reg.all() {
            assert!(!t.label.is_empty());
            assert!(!t.display_name.is_empty());
        }
    }

    #[test]
    fn dock_params_defaults_are_sensible() {
        let p = DockParams::default();
        assert!(p.grid.is_none());
        assert_eq!(p.cluster_rmsd, 2.0);
        assert_eq!(p.n_runs, 8);
    }
}
