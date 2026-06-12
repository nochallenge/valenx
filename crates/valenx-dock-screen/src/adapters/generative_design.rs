//! Feature 27 — generative-protein-design subprocess adapter.
//!
//! **Adapter only — no neural network is reimplemented.** This module
//! wraps the external generative-protein-design tools: RFdiffusion
//! (backbone generation), ProteinMPNN and ESM-IF (inverse folding —
//! sequence design for a given backbone), and Chroma (a generative
//! protein model). All are trained deep-learning models; Valenx never
//! reimplements them — see the [`crate::adapters`] module docs.
//!
//! [`run_generative_design`] takes a typed
//! [`GenerativeDesignRequest`], locates the chosen tool on `PATH`, and
//! either returns the [`crate::adapters::AdapterCommand`]
//! that would run it or a
//! [`crate::error::DockScreenError::ToolNotAvailable`].

use std::path::PathBuf;

use crate::adapters::common::{find_executable, AdapterCommand, ToolStatus};
use crate::error::{DockScreenError, Result};

/// An external generative-protein-design tool.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GenerativeTool {
    /// RFdiffusion — diffusion-based protein backbone generation.
    RfDiffusion,
    /// ProteinMPNN — message-passing inverse folding.
    ProteinMpnn,
    /// ESM-IF — ESM inverse folding.
    EsmIf,
    /// Chroma — a programmable generative protein model.
    Chroma,
}

impl GenerativeTool {
    /// A short stable label.
    pub fn label(self) -> &'static str {
        match self {
            GenerativeTool::RfDiffusion => "rfdiffusion",
            GenerativeTool::ProteinMpnn => "proteinmpnn",
            GenerativeTool::EsmIf => "esm_if",
            GenerativeTool::Chroma => "chroma",
        }
    }

    /// A human-readable display name.
    pub fn display_name(self) -> &'static str {
        match self {
            GenerativeTool::RfDiffusion => "RFdiffusion",
            GenerativeTool::ProteinMpnn => "ProteinMPNN",
            GenerativeTool::EsmIf => "ESM-IF",
            GenerativeTool::Chroma => "Chroma",
        }
    }

    /// `true` if the tool *generates a backbone* (a 3D structure);
    /// `false` if it *designs a sequence* for an existing backbone
    /// (inverse folding).
    pub fn generates_backbone(self) -> bool {
        match self {
            GenerativeTool::RfDiffusion | GenerativeTool::Chroma => true,
            GenerativeTool::ProteinMpnn | GenerativeTool::EsmIf => false,
        }
    }

    /// Candidate binary names, in preference order.
    fn binary_candidates(self) -> &'static [&'static str] {
        match self {
            GenerativeTool::RfDiffusion => &["run_inference.py", "rfdiffusion"],
            GenerativeTool::ProteinMpnn => &["protein_mpnn_run.py", "proteinmpnn"],
            GenerativeTool::EsmIf => &["esm-if", "esm_inverse_folding"],
            GenerativeTool::Chroma => &["chroma"],
        }
    }

    /// An install hint surfaced when the tool is missing.
    fn install_hint(self) -> &'static str {
        match self {
            GenerativeTool::RfDiffusion => {
                "install RFdiffusion from github.com/RosettaCommons/RFdiffusion"
            }
            GenerativeTool::ProteinMpnn => {
                "install ProteinMPNN from github.com/dauparas/ProteinMPNN"
            }
            GenerativeTool::EsmIf => "install ESM-IF (`pip install fair-esm`)",
            GenerativeTool::Chroma => "install Chroma (`pip install generate-chroma`)",
        }
    }
}

/// A typed generative-protein-design request.
#[derive(Clone, Debug, PartialEq)]
pub struct GenerativeDesignRequest {
    /// Which external tool to run.
    pub tool: GenerativeTool,
    /// Input structure (a PDB) — the scaffold / motif for a
    /// backbone-generating tool, or the backbone to redesign for an
    /// inverse-folding tool.
    pub input_structure: PathBuf,
    /// Directory the tool should write its designs into.
    pub output_dir: PathBuf,
    /// Number of designs (backbones or sequences) to generate.
    pub num_designs: u32,
    /// Sampling temperature (where the tool supports it).
    pub temperature: f64,
}

impl GenerativeDesignRequest {
    /// A request with default settings (8 designs, temperature 0.1).
    pub fn new(
        tool: GenerativeTool,
        input_structure: impl Into<PathBuf>,
        output_dir: impl Into<PathBuf>,
    ) -> Self {
        GenerativeDesignRequest {
            tool,
            input_structure: input_structure.into(),
            output_dir: output_dir.into(),
            num_designs: 8,
            temperature: 0.1,
        }
    }
}

/// The result of preparing a generative-design job.
#[derive(Clone, Debug, PartialEq)]
pub struct GenerativeDesignResult {
    /// The tool that was selected.
    pub tool: GenerativeTool,
    /// The subprocess command to run (the caller executes it).
    pub command: AdapterCommand,
    /// The directory the designs will be written to.
    pub output_dir: PathBuf,
}

/// Feature 27 — prepare a generative-protein-design subprocess job.
///
/// Locates the requested tool on `PATH`; returns the [`AdapterCommand`]
/// when present, or [`DockScreenError::ToolNotAvailable`] when absent.
///
/// **No neural network is run here, and none is reimplemented.**
pub fn run_generative_design(request: &GenerativeDesignRequest) -> Result<GenerativeDesignResult> {
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
            "{} generative design from {}",
            request.tool.display_name(),
            request.input_structure.display()
        ),
    );
    Ok(GenerativeDesignResult {
        tool: request.tool,
        command,
        output_dir: request.output_dir.clone(),
    })
}

/// Build the per-tool argument vector.
fn build_args(request: &GenerativeDesignRequest) -> Vec<String> {
    let input = request.input_structure.display().to_string();
    let out = request.output_dir.display().to_string();
    match request.tool {
        GenerativeTool::RfDiffusion => vec![
            format!("inference.input_pdb={input}"),
            format!("inference.output_prefix={out}/design"),
            format!("inference.num_designs={}", request.num_designs),
        ],
        GenerativeTool::ProteinMpnn => vec![
            "--pdb_path".into(),
            input,
            "--out_folder".into(),
            out,
            "--num_seq_per_target".into(),
            request.num_designs.to_string(),
            "--sampling_temp".into(),
            format!("{:.3}", request.temperature),
        ],
        GenerativeTool::EsmIf => vec![
            "--pdb".into(),
            input,
            "--outpath".into(),
            out,
            "--num-samples".into(),
            request.num_designs.to_string(),
            "--temperature".into(),
            format!("{:.3}", request.temperature),
        ],
        GenerativeTool::Chroma => vec![
            "sample".into(),
            "--init".into(),
            input,
            "--out".into(),
            out,
            "--samples".into(),
            request.num_designs.to_string(),
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_tool_returns_tool_not_available() {
        for tool in [
            GenerativeTool::RfDiffusion,
            GenerativeTool::ProteinMpnn,
            GenerativeTool::EsmIf,
            GenerativeTool::Chroma,
        ] {
            let req = GenerativeDesignRequest::new(tool, "scaffold.pdb", "out");
            match run_generative_design(&req) {
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
    fn backbone_vs_sequence_classification() {
        assert!(GenerativeTool::RfDiffusion.generates_backbone());
        assert!(GenerativeTool::Chroma.generates_backbone());
        // Inverse-folding tools design a sequence, not a backbone.
        assert!(!GenerativeTool::ProteinMpnn.generates_backbone());
        assert!(!GenerativeTool::EsmIf.generates_backbone());
    }

    #[test]
    fn labels_are_stable() {
        assert_eq!(GenerativeTool::RfDiffusion.label(), "rfdiffusion");
        assert_eq!(GenerativeTool::ProteinMpnn.display_name(), "ProteinMPNN");
    }

    #[test]
    fn request_defaults() {
        let req = GenerativeDesignRequest::new(GenerativeTool::ProteinMpnn, "in.pdb", "out");
        assert_eq!(req.num_designs, 8);
        assert!((req.temperature - 0.1).abs() < 1e-12);
    }

    #[test]
    fn build_args_embeds_paths_and_counts() {
        let req =
            GenerativeDesignRequest::new(GenerativeTool::ProteinMpnn, "backbone.pdb", "designs");
        let args = build_args(&req);
        assert!(args.iter().any(|a| a.contains("backbone.pdb")));
        assert!(args.iter().any(|a| a.contains("designs")));
        // The design count must reach the command line.
        assert!(args.iter().any(|a| a == "8"));
    }

    #[test]
    fn rfdiffusion_args_use_its_hydra_style_keys() {
        let req = GenerativeDesignRequest::new(GenerativeTool::RfDiffusion, "motif.pdb", "out");
        let args = build_args(&req);
        // RFdiffusion uses `key=value` Hydra-style arguments.
        assert!(args.iter().any(|a| a.starts_with("inference.input_pdb=")));
        assert!(args.iter().any(|a| a.starts_with("inference.num_designs=")));
    }
}
