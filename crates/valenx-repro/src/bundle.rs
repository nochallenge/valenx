//! The reproducibility bundle and its artifacts.

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::error::ReproError;

/// Whether an [`Artifact`] is an input to, or an output of, the study.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArtifactRole {
    /// An input (e.g. a reference sequence, a parameter file).
    Input,
    /// An output / result (e.g. a predicted structure, a report).
    Output,
}

impl ArtifactRole {
    fn tag(self) -> &'static str {
        match self {
            ArtifactRole::Input => "input",
            ArtifactRole::Output => "output",
        }
    }
}

/// A single data artifact, identified by name and content hash.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Artifact {
    /// File name or logical identifier.
    pub name: String,
    /// Input or output.
    pub role: ArtifactRole,
    /// Lower-case hex SHA-256 of the artifact's bytes.
    pub sha256: String,
    /// Size in bytes.
    pub size_bytes: u64,
}

impl Artifact {
    /// Build an artifact by hashing its bytes with SHA-256.
    pub fn from_bytes(name: impl Into<String>, role: ArtifactRole, bytes: &[u8]) -> Self {
        let mut hasher = Sha256::new();
        hasher.update(bytes);
        Self {
            name: name.into(),
            role,
            sha256: hex(&hasher.finalize()),
            size_bytes: bytes.len() as u64,
        }
    }
}

/// A named parameter recorded for reproducibility.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Parameter {
    /// Parameter name.
    pub name: String,
    /// Parameter value, rendered as text.
    pub value: String,
}

impl Parameter {
    /// A new parameter.
    pub fn new(name: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            value: value.into(),
        }
    }
}

/// One step of the workflow that produced the results.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProvenanceStep {
    /// Execution order; must be unique within a bundle.
    pub ordinal: u32,
    /// The tool or function invoked.
    pub tool: String,
    /// The tool's version string.
    pub version: String,
    /// The command line / arguments used.
    pub args: String,
}

impl ProvenanceStep {
    /// A new provenance step.
    pub fn new(
        ordinal: u32,
        tool: impl Into<String>,
        version: impl Into<String>,
        args: impl Into<String>,
    ) -> Self {
        Self {
            ordinal,
            tool: tool.into(),
            version: version.into(),
            args: args.into(),
        }
    }
}

/// A piece of software, with its version, in the bundle's manifest.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SoftwareRef {
    /// Software name.
    pub name: String,
    /// Version string.
    pub version: String,
}

impl SoftwareRef {
    /// A new software reference.
    pub fn new(name: impl Into<String>, version: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            version: version.into(),
        }
    }
}

/// A complete, verifiable reproducibility bundle for a computational study.
///
/// Build it incrementally, then read off [`ReproBundle::fingerprint`] (a
/// deterministic SHA-256 over the canonical contents — any change to any
/// artifact, parameter, software entry or step changes it) and the templated
/// [`crate::methods_scaffold`] / [`crate::abstract_scaffold`].
///
/// This type is **export-only**: it writes no files and performs no network
/// action. It is `serde`-serializable; persist it yourself.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReproBundle {
    /// Study title.
    pub title: String,
    /// Short description / abstract seed.
    pub description: String,
    /// Input artifacts.
    pub inputs: Vec<Artifact>,
    /// Recorded parameters.
    pub parameters: Vec<Parameter>,
    /// Ordered workflow provenance.
    pub steps: Vec<ProvenanceStep>,
    /// Output / result artifacts.
    pub outputs: Vec<Artifact>,
    /// Software / version manifest.
    pub software: Vec<SoftwareRef>,
}

impl ReproBundle {
    /// Start a new bundle. `title` and `description` must be non-empty.
    pub fn new(
        title: impl Into<String>,
        description: impl Into<String>,
    ) -> Result<Self, ReproError> {
        let title = title.into();
        let description = description.into();
        if title.trim().is_empty() {
            return Err(ReproError::Empty { field: "title" });
        }
        if description.trim().is_empty() {
            return Err(ReproError::Empty {
                field: "description",
            });
        }
        Ok(Self {
            title,
            description,
            inputs: Vec::new(),
            parameters: Vec::new(),
            steps: Vec::new(),
            outputs: Vec::new(),
            software: Vec::new(),
        })
    }

    /// Add an artifact, filed under inputs or outputs by its `role`.
    pub fn with_artifact(mut self, artifact: Artifact) -> Self {
        match artifact.role {
            ArtifactRole::Input => self.inputs.push(artifact),
            ArtifactRole::Output => self.outputs.push(artifact),
        }
        self
    }

    /// Add a parameter.
    pub fn with_parameter(mut self, parameter: Parameter) -> Self {
        self.parameters.push(parameter);
        self
    }

    /// Add a software reference.
    pub fn with_software(mut self, software: SoftwareRef) -> Self {
        self.software.push(software);
        self
    }

    /// Add a provenance step, rejecting a duplicate ordinal.
    pub fn with_step(mut self, step: ProvenanceStep) -> Result<Self, ReproError> {
        if self.steps.iter().any(|s| s.ordinal == step.ordinal) {
            return Err(ReproError::DuplicateStepOrdinal(step.ordinal));
        }
        self.steps.push(step);
        Ok(self)
    }

    /// The deterministic SHA-256 fingerprint of the whole bundle (lower-case
    /// hex). Canonicalised by sorting artifacts / parameters / software by
    /// name (so insertion order does not matter) and steps by ordinal. Any
    /// content change flips it — this is the tamper-evidence root.
    pub fn fingerprint(&self) -> String {
        let mut buf = String::new();
        buf.push_str("title\0");
        buf.push_str(&self.title);
        buf.push_str("\0description\0");
        buf.push_str(&self.description);

        let mut inputs: Vec<&Artifact> = self.inputs.iter().collect();
        inputs.sort_by(|a, b| a.name.cmp(&b.name));
        let mut outputs: Vec<&Artifact> = self.outputs.iter().collect();
        outputs.sort_by(|a, b| a.name.cmp(&b.name));
        for a in inputs.into_iter().chain(outputs) {
            push_artifact(&mut buf, a);
        }

        let mut params: Vec<&Parameter> = self.parameters.iter().collect();
        params.sort_by(|a, b| a.name.cmp(&b.name));
        for p in params {
            buf.push_str("\0param\0");
            buf.push_str(&p.name);
            buf.push('\0');
            buf.push_str(&p.value);
        }

        let mut sw: Vec<&SoftwareRef> = self.software.iter().collect();
        sw.sort_by(|a, b| a.name.cmp(&b.name));
        for s in sw {
            buf.push_str("\0sw\0");
            buf.push_str(&s.name);
            buf.push('\0');
            buf.push_str(&s.version);
        }

        let mut steps: Vec<&ProvenanceStep> = self.steps.iter().collect();
        steps.sort_by_key(|s| s.ordinal);
        for st in steps {
            buf.push_str("\0step\0");
            let _ = write!(buf, "{}", st.ordinal);
            buf.push('\0');
            buf.push_str(&st.tool);
            buf.push('\0');
            buf.push_str(&st.version);
            buf.push('\0');
            buf.push_str(&st.args);
        }

        let mut hasher = Sha256::new();
        hasher.update(buf.as_bytes());
        hex(&hasher.finalize())
    }

    /// A short human-readable manifest summary (counts + the fingerprint).
    pub fn manifest_summary(&self) -> String {
        format!(
            "Reproducibility bundle: {title}\n  inputs: {ni}  outputs: {no}  parameters: {np}  steps: {ns}  software: {nsw}\n  fingerprint (sha256): {fp}",
            title = self.title,
            ni = self.inputs.len(),
            no = self.outputs.len(),
            np = self.parameters.len(),
            ns = self.steps.len(),
            nsw = self.software.len(),
            fp = self.fingerprint(),
        )
    }
}

fn push_artifact(buf: &mut String, a: &Artifact) {
    buf.push_str("\0artifact\0");
    buf.push_str(a.role.tag());
    buf.push('\0');
    buf.push_str(&a.name);
    buf.push('\0');
    buf.push_str(&a.sha256);
    buf.push('\0');
    let _ = write!(buf, "{}", a.size_bytes);
}

fn hex(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        let _ = write!(s, "{b:02x}");
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> ReproBundle {
        ReproBundle::new("MSTN variant study", "Effect of a nonsense MSTN variant")
            .unwrap()
            .with_artifact(Artifact::from_bytes(
                "mstn_cds.fasta",
                ArtifactRole::Input,
                b"ATGCAA",
            ))
            .with_parameter(Parameter::new("genetic_code", "1"))
            .with_software(SoftwareRef::new("valenx-bioseq", "0.1.0"))
            .with_step(ProvenanceStep::new(1, "translate", "0.1.0", "frame +1"))
            .unwrap()
            .with_artifact(Artifact::from_bytes(
                "protein.txt",
                ArtifactRole::Output,
                b"MQ",
            ))
    }

    #[test]
    fn artifact_hashes_known_value() {
        // SHA-256("abc") is a well-known constant.
        let a = Artifact::from_bytes("x.txt", ArtifactRole::Input, b"abc");
        assert_eq!(
            a.sha256,
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(a.size_bytes, 3);
    }

    #[test]
    fn fingerprint_is_deterministic_and_64_hex() {
        assert_eq!(sample().fingerprint(), sample().fingerprint());
        assert_eq!(sample().fingerprint().len(), 64);
        assert!(sample()
            .fingerprint()
            .chars()
            .all(|c| c.is_ascii_hexdigit()));
    }

    #[test]
    fn fingerprint_is_insertion_order_independent() {
        let a = ReproBundle::new("t", "d")
            .unwrap()
            .with_artifact(Artifact::from_bytes("a", ArtifactRole::Input, b"1"))
            .with_artifact(Artifact::from_bytes("b", ArtifactRole::Input, b"2"));
        let b = ReproBundle::new("t", "d")
            .unwrap()
            .with_artifact(Artifact::from_bytes("b", ArtifactRole::Input, b"2"))
            .with_artifact(Artifact::from_bytes("a", ArtifactRole::Input, b"1"));
        assert_eq!(a.fingerprint(), b.fingerprint());
    }

    #[test]
    fn fingerprint_changes_when_anything_changes() {
        let base = sample().fingerprint();
        assert_ne!(
            base,
            sample()
                .with_parameter(Parameter::new("extra", "x"))
                .fingerprint()
        );
    }

    #[test]
    fn rejects_empty_title_and_duplicate_ordinal() {
        assert!(ReproBundle::new("   ", "d").is_err());
        assert!(ReproBundle::new("t", "").is_err());
        let r = ReproBundle::new("t", "d")
            .unwrap()
            .with_step(ProvenanceStep::new(1, "a", "1", ""))
            .unwrap()
            .with_step(ProvenanceStep::new(1, "b", "1", ""));
        assert_eq!(r.unwrap_err(), ReproError::DuplicateStepOrdinal(1));
    }

    #[test]
    fn manifest_summary_mentions_counts_and_fingerprint() {
        let s = sample().manifest_summary();
        assert!(s.contains(&sample().fingerprint()));
        assert!(s.contains("inputs: 1"));
        assert!(s.contains("outputs: 1"));
    }
}
