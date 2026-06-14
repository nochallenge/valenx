//! Templated methods / abstract scaffolds generated from a bundle's
//! provenance.
//!
//! No language model is involved — this is deterministic fill-in-the-blanks
//! templating from the recorded facts, with `[...]` placeholders where the
//! author must supply narrative.

use std::fmt::Write as _;

use crate::bundle::{ProvenanceStep, ReproBundle};

/// A templated **Methods** section built from the bundle's software,
/// inputs, parameters and ordered steps, ending with the reproducibility
/// fingerprint.
pub fn methods_scaffold(bundle: &ReproBundle) -> String {
    let mut out = String::from("## Methods\n\n");

    if bundle.software.is_empty() {
        out.push_str("Software: [list the tools and versions used].\n\n");
    } else {
        let mut sw: Vec<String> = bundle
            .software
            .iter()
            .map(|s| format!("{} {}", s.name, s.version))
            .collect();
        sw.sort();
        out.push_str("All analyses were performed with ");
        out.push_str(&sw.join(", "));
        out.push_str(".\n\n");
    }

    if !bundle.inputs.is_empty() {
        let names: Vec<&str> = bundle.inputs.iter().map(|a| a.name.as_str()).collect();
        out.push_str("Inputs: ");
        out.push_str(&names.join(", "));
        out.push_str(".\n\n");
    }

    if !bundle.parameters.is_empty() {
        let params: Vec<String> = bundle
            .parameters
            .iter()
            .map(|p| format!("{} = {}", p.name, p.value))
            .collect();
        out.push_str("Parameters: ");
        out.push_str(&params.join("; "));
        out.push_str(".\n\n");
    }

    if !bundle.steps.is_empty() {
        out.push_str("Workflow:\n");
        let mut steps: Vec<&ProvenanceStep> = bundle.steps.iter().collect();
        steps.sort_by_key(|s| s.ordinal);
        for s in steps {
            let _ = writeln!(
                out,
                "  {}. {} ({}): {}",
                s.ordinal, s.tool, s.version, s.args
            );
        }
        out.push('\n');
    }

    let _ = writeln!(
        out,
        "Reproducibility: bundle fingerprint (SHA-256) {}.",
        bundle.fingerprint()
    );
    out
}

/// A templated **Abstract** skeleton: recorded facts where available, `[...]`
/// placeholders for the narrative — and an explicit reminder that an in-silico
/// result is a hypothesis, not a validated finding.
pub fn abstract_scaffold(bundle: &ReproBundle) -> String {
    format!(
        "## Abstract\n\n\
         **Background.** [State the question and why it matters.]\n\n\
         **Methods.** {desc} Analyses used {nsw} software tool(s) over {ns} workflow step(s); see Methods.\n\n\
         **Results.** [State the key in-silico finding. Note: this is a computational hypothesis, not a validated result.]\n\n\
         **Conclusion.** [State what should be tested experimentally to validate it.]\n",
        desc = bundle.description,
        nsw = bundle.software.len(),
        ns = bundle.steps.len(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bundle::{Artifact, ArtifactRole, Parameter, ProvenanceStep, SoftwareRef};

    fn bundle() -> ReproBundle {
        ReproBundle::new("Study", "A test study.")
            .unwrap()
            .with_software(SoftwareRef::new("valenx-bioseq", "0.1.0"))
            .with_parameter(Parameter::new("code", "1"))
            .with_artifact(Artifact::from_bytes("in.fa", ArtifactRole::Input, b"ATG"))
            .with_step(ProvenanceStep::new(1, "translate", "0.1.0", "frame +1"))
            .unwrap()
    }

    #[test]
    fn methods_lists_software_steps_and_fingerprint() {
        let m = methods_scaffold(&bundle());
        assert!(m.contains("valenx-bioseq 0.1.0"));
        assert!(m.contains("translate"));
        assert!(m.contains("in.fa"));
        assert!(m.contains(&bundle().fingerprint()));
    }

    #[test]
    fn abstract_has_sections_and_honest_caveat() {
        let a = abstract_scaffold(&bundle());
        for section in ["Background", "Methods", "Results", "Conclusion"] {
            assert!(a.contains(section), "missing {section}");
        }
        // The honesty line must survive: a hypothesis, not a validated result.
        assert!(a.contains("hypothesis"));
    }
}
