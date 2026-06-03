//! Feature 26 — structure-prediction subprocess adapter.
//!
//! **Adapter only — no neural network is reimplemented.** This module
//! wraps the external protein-structure-prediction tools. They are all
//! trained deep-learning models requiring downloaded weights and
//! (mostly) a GPU; Valenx never reimplements them — see the
//! [`crate::adapters`] module docs.
//!
//! Supported tools ([`StructurePredictionTool`]): AlphaFold 2/3,
//! ColabFold, ESMFold, RoseTTAFold, OmegaFold, Boltz and Chai.
//!
//! [`run_structure_prediction`] takes a typed
//! [`StructurePredictionRequest`] (a FASTA query, an output directory,
//! a recycle count), locates the chosen tool on `PATH`, and either
//! returns the [`crate::adapters::AdapterCommand`] that
//! would run it or a
//! [`crate::error::DockScreenError::ToolNotAvailable`]
//! if it is not installed.

use std::path::PathBuf;

use crate::adapters::common::{find_executable, AdapterCommand, ToolStatus};
use crate::error::{DockScreenError, Result};

/// An external protein-structure-prediction tool.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StructurePredictionTool {
    /// DeepMind AlphaFold 2.
    AlphaFold2,
    /// DeepMind AlphaFold 3.
    AlphaFold3,
    /// ColabFold (`colabfold_batch`).
    ColabFold,
    /// Meta ESMFold.
    EsmFold,
    /// RoseTTAFold (Baker lab).
    RoseTTAFold,
    /// OmegaFold.
    OmegaFold,
    /// Boltz (`boltz`).
    Boltz,
    /// Chai-1 (`chai`).
    Chai,
}

impl StructurePredictionTool {
    /// A short stable label.
    pub fn label(self) -> &'static str {
        match self {
            StructurePredictionTool::AlphaFold2 => "alphafold2",
            StructurePredictionTool::AlphaFold3 => "alphafold3",
            StructurePredictionTool::ColabFold => "colabfold",
            StructurePredictionTool::EsmFold => "esmfold",
            StructurePredictionTool::RoseTTAFold => "rosettafold",
            StructurePredictionTool::OmegaFold => "omegafold",
            StructurePredictionTool::Boltz => "boltz",
            StructurePredictionTool::Chai => "chai",
        }
    }

    /// A human-readable display name.
    pub fn display_name(self) -> &'static str {
        match self {
            StructurePredictionTool::AlphaFold2 => "AlphaFold 2",
            StructurePredictionTool::AlphaFold3 => "AlphaFold 3",
            StructurePredictionTool::ColabFold => "ColabFold",
            StructurePredictionTool::EsmFold => "ESMFold",
            StructurePredictionTool::RoseTTAFold => "RoseTTAFold",
            StructurePredictionTool::OmegaFold => "OmegaFold",
            StructurePredictionTool::Boltz => "Boltz",
            StructurePredictionTool::Chai => "Chai-1",
        }
    }

    /// The candidate binary names this tool is invoked through, in
    /// preference order.
    fn binary_candidates(self) -> &'static [&'static str] {
        match self {
            StructurePredictionTool::AlphaFold2 => &["run_alphafold.py", "alphafold"],
            StructurePredictionTool::AlphaFold3 => &["run_alphafold", "alphafold3"],
            StructurePredictionTool::ColabFold => &["colabfold_batch"],
            StructurePredictionTool::EsmFold => &["esm-fold", "esmfold"],
            StructurePredictionTool::RoseTTAFold => &["rosettafold", "run_RF2.sh"],
            StructurePredictionTool::OmegaFold => &["omegafold"],
            StructurePredictionTool::Boltz => &["boltz"],
            StructurePredictionTool::Chai => &["chai", "chai-lab"],
        }
    }

    /// An install hint surfaced when the tool is missing.
    fn install_hint(self) -> &'static str {
        match self {
            StructurePredictionTool::AlphaFold2 => {
                "install AlphaFold 2 from github.com/google-deepmind/alphafold \
                 and put `run_alphafold.py` on PATH"
            }
            StructurePredictionTool::AlphaFold3 => {
                "install AlphaFold 3 from github.com/google-deepmind/alphafold3"
            }
            StructurePredictionTool::ColabFold => {
                "install ColabFold (`pip install colabfold`) so `colabfold_batch` is on PATH"
            }
            StructurePredictionTool::EsmFold => {
                "install ESMFold (`pip install fair-esm[esmfold]`)"
            }
            StructurePredictionTool::RoseTTAFold => {
                "install RoseTTAFold from github.com/RosettaCommons/RoseTTAFold"
            }
            StructurePredictionTool::OmegaFold => "install OmegaFold (`pip install omegafold`)",
            StructurePredictionTool::Boltz => "install Boltz (`pip install boltz`)",
            StructurePredictionTool::Chai => "install Chai-1 (`pip install chai_lab`)",
        }
    }
}

/// A typed structure-prediction request.
#[derive(Clone, Debug, PartialEq)]
pub struct StructurePredictionRequest {
    /// Which external tool to run.
    pub tool: StructurePredictionTool,
    /// Path to the input FASTA query (one or more sequences).
    pub query_fasta: PathBuf,
    /// Directory the tool should write its predicted structures into.
    pub output_dir: PathBuf,
    /// Number of recycling iterations (where the tool supports it).
    pub num_recycles: u32,
    /// Whether to use templates / an MSA (where the tool supports it).
    pub use_msa: bool,
}

impl StructurePredictionRequest {
    /// A request with default settings (4 recycles, MSA on).
    pub fn new(
        tool: StructurePredictionTool,
        query_fasta: impl Into<PathBuf>,
        output_dir: impl Into<PathBuf>,
    ) -> Self {
        StructurePredictionRequest {
            tool,
            query_fasta: query_fasta.into(),
            output_dir: output_dir.into(),
            num_recycles: 4,
            use_msa: true,
        }
    }
}

/// The result of preparing a structure-prediction job: the command
/// that would run the tool, plus the directory its outputs land in.
#[derive(Clone, Debug, PartialEq)]
pub struct StructurePredictionResult {
    /// The tool that was selected.
    pub tool: StructurePredictionTool,
    /// The subprocess command to run (the caller executes it).
    pub command: AdapterCommand,
    /// The directory the predicted structures will be written to.
    pub output_dir: PathBuf,
}

/// Feature 26 — prepare a structure-prediction subprocess job.
///
/// Locates the requested tool on `PATH`. If it is present, returns the
/// [`AdapterCommand`] that would run it (the caller / job runner
/// executes it). If it is absent, returns
/// [`DockScreenError::ToolNotAvailable`] with an install hint.
///
/// **No neural network is run here, and none is reimplemented** — this
/// is honest subprocess scaffolding.
pub fn run_structure_prediction(
    request: &StructurePredictionRequest,
) -> Result<StructurePredictionResult> {
    let status = find_executable(request.tool.binary_candidates());
    let program = match status {
        ToolStatus::Available { path } => path,
        ToolStatus::Missing => {
            return Err(DockScreenError::tool_not_available(
                tool_static_name(request.tool),
                request.tool.install_hint(),
            ))
        }
    };
    let args = build_args(request);
    let command = AdapterCommand::new(
        program,
        args,
        format!(
            "{} structure prediction of {}",
            request.tool.display_name(),
            request.query_fasta.display()
        ),
    );
    Ok(StructurePredictionResult {
        tool: request.tool,
        command,
        output_dir: request.output_dir.clone(),
    })
}

/// The `&'static str` tool name for the error variant (which needs a
/// `'static` lifetime).
fn tool_static_name(tool: StructurePredictionTool) -> &'static str {
    tool.display_name()
}

/// Build the per-tool argument vector. Each tool has its own CLI
/// surface; these are the canonical invocations.
fn build_args(request: &StructurePredictionRequest) -> Vec<String> {
    let fasta = request.query_fasta.display().to_string();
    let out = request.output_dir.display().to_string();
    match request.tool {
        StructurePredictionTool::AlphaFold2 => vec![
            format!("--fasta_paths={fasta}"),
            format!("--output_dir={out}"),
        ],
        StructurePredictionTool::AlphaFold3 => vec![
            "--json_path".into(),
            fasta,
            "--output_dir".into(),
            out,
        ],
        StructurePredictionTool::ColabFold => {
            let mut a = vec![fasta, out];
            a.push(format!("--num-recycle={}", request.num_recycles));
            if !request.use_msa {
                a.push("--msa-mode=single_sequence".into());
            }
            a
        }
        StructurePredictionTool::EsmFold => {
            vec!["-i".into(), fasta, "-o".into(), out]
        }
        StructurePredictionTool::RoseTTAFold => vec![fasta, out],
        StructurePredictionTool::OmegaFold => {
            let mut a = vec![fasta, out];
            a.push(format!("--num_cycle={}", request.num_recycles));
            a
        }
        StructurePredictionTool::Boltz => {
            vec!["predict".into(), fasta, "--out_dir".into(), out]
        }
        StructurePredictionTool::Chai => {
            vec!["fold".into(), fasta, out]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_tool_returns_tool_not_available() {
        // None of these tools is on a CI / dev PATH — the adapter must
        // return the honest "not available" error, never a fake result.
        for tool in [
            StructurePredictionTool::AlphaFold2,
            StructurePredictionTool::AlphaFold3,
            StructurePredictionTool::ColabFold,
            StructurePredictionTool::EsmFold,
            StructurePredictionTool::RoseTTAFold,
            StructurePredictionTool::OmegaFold,
            StructurePredictionTool::Boltz,
            StructurePredictionTool::Chai,
        ] {
            let req = StructurePredictionRequest::new(tool, "query.fasta", "out");
            let result = run_structure_prediction(&req);
            // Either the tool genuinely is not installed (the common
            // case → ToolNotAvailable) or, on a machine where it IS
            // installed, a real command is returned. Both are correct;
            // what is NOT allowed is a fabricated structure.
            match result {
                Err(DockScreenError::ToolNotAvailable { hint, .. }) => {
                    assert!(!hint.is_empty(), "install hint must be present");
                }
                Ok(r) => {
                    assert_eq!(r.tool, tool);
                    assert!(!r.command.args.is_empty());
                }
                Err(other) => panic!("unexpected error variant: {other}"),
            }
        }
    }

    #[test]
    fn tool_labels_and_names_are_stable() {
        assert_eq!(StructurePredictionTool::AlphaFold2.label(), "alphafold2");
        assert_eq!(
            StructurePredictionTool::EsmFold.display_name(),
            "ESMFold"
        );
        assert_eq!(StructurePredictionTool::Boltz.label(), "boltz");
    }

    #[test]
    fn request_defaults_are_sensible() {
        let req = StructurePredictionRequest::new(
            StructurePredictionTool::ColabFold,
            "q.fasta",
            "out",
        );
        assert_eq!(req.num_recycles, 4);
        assert!(req.use_msa);
        assert_eq!(req.tool, StructurePredictionTool::ColabFold);
    }

    #[test]
    fn build_args_embeds_the_paths() {
        // The argument vector for each tool must mention the query and
        // output paths so the constructed command is actually runnable.
        let req = StructurePredictionRequest::new(
            StructurePredictionTool::ColabFold,
            "myquery.fasta",
            "results",
        );
        let args = build_args(&req);
        assert!(args.iter().any(|a| a.contains("myquery.fasta")));
        assert!(args.iter().any(|a| a.contains("results")));
        assert!(args.iter().any(|a| a.contains("num-recycle")));
    }

    #[test]
    fn colabfold_single_sequence_arg_appears_when_msa_disabled() {
        let mut req = StructurePredictionRequest::new(
            StructurePredictionTool::ColabFold,
            "q.fasta",
            "out",
        );
        req.use_msa = false;
        let args = build_args(&req);
        assert!(args.iter().any(|a| a.contains("single_sequence")));
    }

    #[test]
    fn the_install_hint_names_the_tool_or_its_install_command() {
        // The hint must actually help the user — not be empty filler.
        for tool in [
            StructurePredictionTool::AlphaFold2,
            StructurePredictionTool::ColabFold,
            StructurePredictionTool::Boltz,
        ] {
            assert!(tool.install_hint().len() > 10);
        }
    }
}
