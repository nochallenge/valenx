//! Feature 28 — neural-network docking subprocess adapter.
//!
//! **Adapter only — no neural network is reimplemented.** This module
//! wraps the external deep-learning docking tools: DiffDock (a
//! diffusion-generative docking model) and GNINA (a CNN-rescored
//! AutoDock fork). Both need trained weights and a GPU; Valenx never
//! reimplements them — see the [`crate::adapters`] module docs.
//!
//! Note the contrast with the rest of this crate: classical
//! AutoDock-class docking (the Vina / AutoDock4 scoring functions, the
//! genetic algorithm, Monte Carlo, iterated local search) **is**
//! implemented as a real working v1 in [`crate::score`] and
//! [`crate::search`]. It is only the *neural-network* docking tools
//! that are adapter-only.
//!
//! [`run_nn_docking`] takes a typed [`NnDockingRequest`] (a receptor,
//! a ligand, an output directory), locates the chosen tool on `PATH`,
//! and either returns the
//! [`crate::adapters::AdapterCommand`] that would run
//! it or a
//! [`crate::error::DockScreenError::ToolNotAvailable`].

use std::path::PathBuf;

use crate::adapters::common::{find_executable, AdapterCommand, ToolStatus};
use crate::error::{DockScreenError, Result};

/// An external neural-network docking tool.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NnDockingTool {
    /// DiffDock — diffusion-generative blind docking.
    DiffDock,
    /// GNINA — a CNN-scored AutoDock fork.
    Gnina,
}

impl NnDockingTool {
    /// A short stable label.
    pub fn label(self) -> &'static str {
        match self {
            NnDockingTool::DiffDock => "diffdock",
            NnDockingTool::Gnina => "gnina",
        }
    }

    /// A human-readable display name.
    pub fn display_name(self) -> &'static str {
        match self {
            NnDockingTool::DiffDock => "DiffDock",
            NnDockingTool::Gnina => "GNINA",
        }
    }

    /// Candidate binary names, in preference order.
    fn binary_candidates(self) -> &'static [&'static str] {
        match self {
            NnDockingTool::DiffDock => &["diffdock", "run_diffdock.py"],
            NnDockingTool::Gnina => &["gnina"],
        }
    }

    /// An install hint surfaced when the tool is missing.
    fn install_hint(self) -> &'static str {
        match self {
            NnDockingTool::DiffDock => "install DiffDock from github.com/gcorso/DiffDock",
            NnDockingTool::Gnina => {
                "install GNINA from github.com/gnina/gnina (a CNN-scored AutoDock fork)"
            }
        }
    }
}

/// A typed neural-network docking request.
#[derive(Clone, Debug, PartialEq)]
pub struct NnDockingRequest {
    /// Which external tool to run.
    pub tool: NnDockingTool,
    /// Path to the receptor structure (PDB or PDBQT).
    pub receptor: PathBuf,
    /// Path to the ligand structure (SDF, MOL2 or PDBQT).
    pub ligand: PathBuf,
    /// Directory the tool should write its docked poses into.
    pub output_dir: PathBuf,
    /// Number of poses to generate / keep.
    pub num_poses: u32,
}

impl NnDockingRequest {
    /// A request with default settings (10 poses).
    pub fn new(
        tool: NnDockingTool,
        receptor: impl Into<PathBuf>,
        ligand: impl Into<PathBuf>,
        output_dir: impl Into<PathBuf>,
    ) -> Self {
        NnDockingRequest {
            tool,
            receptor: receptor.into(),
            ligand: ligand.into(),
            output_dir: output_dir.into(),
            num_poses: 10,
        }
    }
}

/// The result of preparing a neural-network docking job.
#[derive(Clone, Debug, PartialEq)]
pub struct NnDockingResult {
    /// The tool that was selected.
    pub tool: NnDockingTool,
    /// The subprocess command to run (the caller executes it).
    pub command: AdapterCommand,
    /// The directory the docked poses will be written to.
    pub output_dir: PathBuf,
}

/// Feature 28 — prepare a neural-network docking subprocess job.
///
/// Locates the requested tool on `PATH`; returns the [`AdapterCommand`]
/// when present, or [`DockScreenError::ToolNotAvailable`] when absent.
///
/// **No neural network is run here, and none is reimplemented.** For
/// real, in-process docking, use the classical AutoDock-class search
/// in [`crate::search`].
pub fn run_nn_docking(request: &NnDockingRequest) -> Result<NnDockingResult> {
    let status = find_executable(request.tool.binary_candidates());
    let program = match status {
        ToolStatus::Available { path } => path,
        ToolStatus::Missing => {
            return Err(DockScreenError::tool_not_available(
                request.tool.display_name(),
                request.tool.install_hint(),
            ))
        }
    };
    let args = build_args(request);
    let command = AdapterCommand::new(
        program,
        args,
        format!(
            "{} docking of {} into {}",
            request.tool.display_name(),
            request.ligand.display(),
            request.receptor.display()
        ),
    );
    Ok(NnDockingResult {
        tool: request.tool,
        command,
        output_dir: request.output_dir.clone(),
    })
}

/// Build the per-tool argument vector.
fn build_args(request: &NnDockingRequest) -> Vec<String> {
    let receptor = request.receptor.display().to_string();
    let ligand = request.ligand.display().to_string();
    let out = request.output_dir.display().to_string();
    match request.tool {
        NnDockingTool::DiffDock => vec![
            "--protein_path".into(),
            receptor,
            "--ligand".into(),
            ligand,
            "--out_dir".into(),
            out,
            "--samples_per_complex".into(),
            request.num_poses.to_string(),
        ],
        NnDockingTool::Gnina => vec![
            "-r".into(),
            receptor,
            "-l".into(),
            ligand,
            "-o".into(),
            format!("{out}/docked.sdf"),
            "--num_modes".into(),
            request.num_poses.to_string(),
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_tool_returns_tool_not_available() {
        for tool in [NnDockingTool::DiffDock, NnDockingTool::Gnina] {
            let req = NnDockingRequest::new(tool, "receptor.pdb", "ligand.sdf", "out");
            match run_nn_docking(&req) {
                Err(DockScreenError::ToolNotAvailable { hint, .. }) => {
                    assert!(!hint.is_empty());
                }
                Ok(r) => {
                    assert_eq!(r.tool, tool);
                    assert!(!r.command.args.is_empty());
                }
                Err(other) => panic!("unexpected error: {other}"),
            }
        }
    }

    #[test]
    fn labels_are_stable() {
        assert_eq!(NnDockingTool::DiffDock.label(), "diffdock");
        assert_eq!(NnDockingTool::Gnina.display_name(), "GNINA");
    }

    #[test]
    fn request_default_pose_count() {
        let req = NnDockingRequest::new(NnDockingTool::DiffDock, "r.pdb", "l.sdf", "o");
        assert_eq!(req.num_poses, 10);
    }

    #[test]
    fn build_args_embeds_all_three_paths() {
        let req = NnDockingRequest::new(
            NnDockingTool::DiffDock,
            "protein.pdb",
            "small_molecule.sdf",
            "results",
        );
        let args = build_args(&req);
        assert!(args.iter().any(|a| a.contains("protein.pdb")));
        assert!(args.iter().any(|a| a.contains("small_molecule.sdf")));
        assert!(args.iter().any(|a| a.contains("results")));
    }

    #[test]
    fn gnina_args_use_short_flags() {
        let req = NnDockingRequest::new(NnDockingTool::Gnina, "r.pdb", "l.sdf", "o");
        let args = build_args(&req);
        // GNINA's CLI is AutoDock-style: -r / -l / -o.
        assert!(args.iter().any(|a| a == "-r"));
        assert!(args.iter().any(|a| a == "-l"));
        assert!(args.iter().any(|a| a == "-o"));
    }

    #[test]
    fn diffdock_args_use_its_long_flags() {
        let req = NnDockingRequest::new(NnDockingTool::DiffDock, "r.pdb", "l.sdf", "o");
        let args = build_args(&req);
        assert!(args.iter().any(|a| a == "--protein_path"));
        assert!(args.iter().any(|a| a == "--samples_per_complex"));
    }
}
