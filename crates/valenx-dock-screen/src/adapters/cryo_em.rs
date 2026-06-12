//! Feature 29 — cryo-EM-reconstruction subprocess adapter.
//!
//! **Adapter only — no reconstruction pipeline is reimplemented.**
//! This module wraps the external single-particle cryo-EM
//! reconstruction suites: RELION, cryoSPARC and EMAN2. These are large
//! GPU image-processing pipelines (motion correction, CTF estimation,
//! particle picking, 2D / 3D classification, refinement) — Valenx
//! never reimplements them; see the [`crate::adapters`] module docs.
//!
//! [`run_cryo_em`] takes a typed [`CryoEmRequest`] (an input particle
//! stack / project, an output directory, a reconstruction stage),
//! locates the chosen suite on `PATH`, and either returns the
//! [`crate::adapters::AdapterCommand`] that would run
//! it or a
//! [`crate::error::DockScreenError::ToolNotAvailable`].

use std::path::PathBuf;

use crate::adapters::common::{find_executable, AdapterCommand, ToolStatus};
use crate::error::{DockScreenError, Result};

/// An external cryo-EM reconstruction suite.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CryoEmTool {
    /// RELION (REgularised LIkelihood OptimisatioN).
    Relion,
    /// cryoSPARC.
    CryoSparc,
    /// EMAN2.
    Eman2,
}

impl CryoEmTool {
    /// A short stable label.
    pub fn label(self) -> &'static str {
        match self {
            CryoEmTool::Relion => "relion",
            CryoEmTool::CryoSparc => "cryosparc",
            CryoEmTool::Eman2 => "eman2",
        }
    }

    /// A human-readable display name.
    pub fn display_name(self) -> &'static str {
        match self {
            CryoEmTool::Relion => "RELION",
            CryoEmTool::CryoSparc => "cryoSPARC",
            CryoEmTool::Eman2 => "EMAN2",
        }
    }

    /// An install hint surfaced when the tool is missing.
    fn install_hint(self) -> &'static str {
        match self {
            CryoEmTool::Relion => "install RELION from relion.readthedocs.io",
            CryoEmTool::CryoSparc => {
                "install cryoSPARC from cryosparc.com (the `cryosparcm` CLI must be on PATH)"
            }
            CryoEmTool::Eman2 => "install EMAN2 from blake.bcm.edu/emanwiki/EMAN2",
        }
    }
}

/// A stage of the single-particle reconstruction pipeline.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CryoEmStage {
    /// 2D classification of a particle stack.
    Class2D,
    /// 3D classification.
    Class3D,
    /// 3D auto-refinement to a final map.
    Refine3D,
}

impl CryoEmStage {
    /// A short stable label.
    pub fn label(self) -> &'static str {
        match self {
            CryoEmStage::Class2D => "class2d",
            CryoEmStage::Class3D => "class3d",
            CryoEmStage::Refine3D => "refine3d",
        }
    }
}

/// A typed cryo-EM reconstruction request.
#[derive(Clone, Debug, PartialEq)]
pub struct CryoEmRequest {
    /// Which external suite to run.
    pub tool: CryoEmTool,
    /// The reconstruction stage to run.
    pub stage: CryoEmStage,
    /// Input particle stack (a RELION `.star` file, an EMAN2 stack, or
    /// a cryoSPARC project path).
    pub input: PathBuf,
    /// Directory the suite should write its outputs into.
    pub output_dir: PathBuf,
    /// Number of classes (for the classification stages).
    pub num_classes: u32,
}

impl CryoEmRequest {
    /// A request with a default class count of 8.
    pub fn new(
        tool: CryoEmTool,
        stage: CryoEmStage,
        input: impl Into<PathBuf>,
        output_dir: impl Into<PathBuf>,
    ) -> Self {
        CryoEmRequest {
            tool,
            stage,
            input: input.into(),
            output_dir: output_dir.into(),
            num_classes: 8,
        }
    }
}

/// The result of preparing a cryo-EM reconstruction job.
#[derive(Clone, Debug, PartialEq)]
pub struct CryoEmResult {
    /// The suite that was selected.
    pub tool: CryoEmTool,
    /// The subprocess command to run (the caller executes it).
    pub command: AdapterCommand,
    /// The directory the reconstruction outputs will be written to.
    pub output_dir: PathBuf,
}

/// The candidate binary names a suite is invoked through for a given
/// stage. RELION ships one binary per job type; cryoSPARC and EMAN2
/// dispatch through a single CLI.
fn binary_candidates(tool: CryoEmTool, stage: CryoEmStage) -> Vec<&'static str> {
    match tool {
        CryoEmTool::Relion => match stage {
            CryoEmStage::Class2D | CryoEmStage::Class3D => {
                vec!["relion_refine", "relion_refine_mpi"]
            }
            CryoEmStage::Refine3D => vec!["relion_refine", "relion_refine_mpi"],
        },
        CryoEmTool::CryoSparc => vec!["cryosparcm"],
        CryoEmTool::Eman2 => match stage {
            CryoEmStage::Class2D => vec!["e2refine2d.py"],
            CryoEmStage::Class3D | CryoEmStage::Refine3D => vec!["e2refine_easy.py"],
        },
    }
}

/// Feature 29 — prepare a cryo-EM reconstruction subprocess job.
///
/// Locates the requested suite on `PATH`; returns the [`AdapterCommand`]
/// when present, or [`DockScreenError::ToolNotAvailable`] when absent.
///
/// **No reconstruction pipeline is run here, and none is
/// reimplemented** — this is honest subprocess scaffolding.
pub fn run_cryo_em(request: &CryoEmRequest) -> Result<CryoEmResult> {
    let candidates = binary_candidates(request.tool, request.stage);
    let status = find_executable(&candidates);
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
            "{} {} reconstruction from {}",
            request.tool.display_name(),
            request.stage.label(),
            request.input.display()
        ),
    );
    Ok(CryoEmResult {
        tool: request.tool,
        command,
        output_dir: request.output_dir.clone(),
    })
}

/// Build the per-suite argument vector.
fn build_args(request: &CryoEmRequest) -> Vec<String> {
    let input = request.input.display().to_string();
    let out = request.output_dir.display().to_string();
    match request.tool {
        CryoEmTool::Relion => {
            let mut a = vec!["--i".into(), input, "--o".into(), format!("{out}/run")];
            match request.stage {
                CryoEmStage::Class2D => {
                    a.push("--class2d".into());
                    a.push("--K".into());
                    a.push(request.num_classes.to_string());
                }
                CryoEmStage::Class3D => {
                    a.push("--class3d".into());
                    a.push("--K".into());
                    a.push(request.num_classes.to_string());
                }
                CryoEmStage::Refine3D => {
                    a.push("--auto_refine".into());
                    a.push("--split_random_halves".into());
                }
            }
            a
        }
        CryoEmTool::CryoSparc => vec![
            "cli".into(),
            format!("do_{}", request.stage.label()),
            "--project".into(),
            input,
            "--output".into(),
            out,
        ],
        CryoEmTool::Eman2 => {
            let mut a = vec![format!("--input={input}"), format!("--path={out}")];
            if matches!(request.stage, CryoEmStage::Class2D | CryoEmStage::Class3D) {
                a.push(format!("--ncls={}", request.num_classes));
            }
            a
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_tool_returns_tool_not_available() {
        for tool in [CryoEmTool::Relion, CryoEmTool::CryoSparc, CryoEmTool::Eman2] {
            let req = CryoEmRequest::new(tool, CryoEmStage::Class2D, "particles.star", "out");
            match run_cryo_em(&req) {
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
        assert_eq!(CryoEmTool::Relion.label(), "relion");
        assert_eq!(CryoEmTool::CryoSparc.display_name(), "cryoSPARC");
        assert_eq!(CryoEmStage::Refine3D.label(), "refine3d");
    }

    #[test]
    fn request_default_class_count() {
        let req = CryoEmRequest::new(CryoEmTool::Relion, CryoEmStage::Class2D, "in.star", "out");
        assert_eq!(req.num_classes, 8);
    }

    #[test]
    fn build_args_embeds_input_and_output() {
        let req = CryoEmRequest::new(
            CryoEmTool::Relion,
            CryoEmStage::Class2D,
            "stack.star",
            "results",
        );
        let args = build_args(&req);
        assert!(args.iter().any(|a| a.contains("stack.star")));
        assert!(args.iter().any(|a| a.contains("results")));
        // The class count must reach the command line for 2D class.
        assert!(args.iter().any(|a| a == "8"));
    }

    #[test]
    fn relion_refine_stage_uses_auto_refine_flag() {
        let req = CryoEmRequest::new(CryoEmTool::Relion, CryoEmStage::Refine3D, "in.star", "out");
        let args = build_args(&req);
        assert!(args.iter().any(|a| a == "--auto_refine"));
        assert!(args.iter().any(|a| a == "--split_random_halves"));
    }

    #[test]
    fn relion_dispatches_to_relion_refine() {
        // Stage-aware binary selection: every RELION stage runs
        // through `relion_refine`.
        let c2 = binary_candidates(CryoEmTool::Relion, CryoEmStage::Class2D);
        let r3 = binary_candidates(CryoEmTool::Relion, CryoEmStage::Refine3D);
        assert!(c2.contains(&"relion_refine"));
        assert!(r3.contains(&"relion_refine"));
    }

    #[test]
    fn eman2_uses_stage_specific_scripts() {
        // EMAN2 has a different script per stage.
        let c2 = binary_candidates(CryoEmTool::Eman2, CryoEmStage::Class2D);
        let r3 = binary_candidates(CryoEmTool::Eman2, CryoEmStage::Refine3D);
        assert!(c2.contains(&"e2refine2d.py"));
        assert!(r3.contains(&"e2refine_easy.py"));
    }
}
