//! SBOL-class genetic-design data model — feature 24.
//!
//! SBOL (the Synthetic Biology Open Language) is the standard
//! interchange format for *genetic designs* — the structure of an
//! engineered DNA construct, as opposed to the kinetics of the cell it
//! runs in. This module provides an SBOL-flavoured data model
//! covering the parts of SBOL that a v1 design tool actually uses.
//!
//! The hierarchy mirrors SBOL's `ComponentDefinition` /
//! `SequenceAnnotation` / `Module` structure:
//!
//! - [`Part`] — an atomic genetic element with a [`PartRole`]
//!   (promoter, RBS, CDS, terminator, …) and an optional DNA
//!   sequence. The SBOL `ComponentDefinition` of a leaf component.
//! - [`SequenceAnnotation`] — a typed feature spanning a sub-range of
//!   a component's sequence.
//! - [`Component`] — an ordered assembly of parts (a transcription
//!   unit, an operon) — a composite `ComponentDefinition` whose
//!   `sequence` is the concatenation of its parts'.
//! - [`Device`] — a named collection of components plus
//!   [`Module`]-style functional interactions, the SBOL `Module`.
//!
//! Sequences use [`valenx_bioseq::Seq`] so the assembly planners and
//! the rest of the Valenx bio stack interoperate directly.

use valenx_bioseq::{Seq, SeqKind};

use crate::error::{Result, SysbioError};

/// The functional role of a genetic [`Part`] — the SBOL / SO role
/// vocabulary, trimmed to the elements a v1 circuit designer uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PartRole {
    /// A promoter — transcription start.
    Promoter,
    /// An operator / transcription-factor binding site.
    Operator,
    /// A ribosome binding site.
    Rbs,
    /// A protein-coding sequence.
    Cds,
    /// A transcription terminator.
    Terminator,
    /// A spacer / scar / assembly junction.
    Spacer,
    /// An origin of replication.
    Origin,
    /// A primer-binding / annotation site.
    PrimerBinding,
    /// Anything else.
    Other,
}

impl PartRole {
    /// The Sequence Ontology-style short tag for this role.
    pub fn so_tag(&self) -> &'static str {
        match self {
            PartRole::Promoter => "SO:promoter",
            PartRole::Operator => "SO:operator",
            PartRole::Rbs => "SO:ribosome_entry_site",
            PartRole::Cds => "SO:CDS",
            PartRole::Terminator => "SO:terminator",
            PartRole::Spacer => "SO:sequence_feature",
            PartRole::Origin => "SO:origin_of_replication",
            PartRole::PrimerBinding => "SO:primer_binding_site",
            PartRole::Other => "SO:sequence_feature",
        }
    }
}

/// An atomic genetic part — the leaf of the design hierarchy.
#[derive(Debug, Clone, PartialEq)]
pub struct Part {
    /// Unique part identifier.
    pub id: String,
    /// Human-readable display name.
    pub name: String,
    /// Functional role.
    pub role: PartRole,
    /// The part's DNA sequence (may be empty for an abstract part).
    pub sequence: Seq,
}

impl Part {
    /// A DNA part from an id, role and IUPAC sequence string.
    pub fn new(id: impl Into<String>, role: PartRole, seq: impl AsRef<[u8]>) -> Result<Self> {
        let id = id.into();
        let sequence = Seq::new(SeqKind::Dna, seq).map_err(|e| {
            SysbioError::invalid_model("sbol", format!("part `{id}` has a bad sequence: {e}"))
        })?;
        Ok(Part {
            name: id.clone(),
            id,
            role,
            sequence,
        })
    }

    /// Builder: set a display name distinct from the id.
    pub fn named(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    /// Length of the part's DNA sequence.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Whether the part has an empty sequence.
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// A typed annotation over a sub-range of a sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct SequenceAnnotation {
    /// Annotation identifier.
    pub id: String,
    /// 0-based inclusive start offset.
    pub start: usize,
    /// 0-based exclusive end offset.
    pub end: usize,
    /// The role of the annotated feature.
    pub role: PartRole,
}

impl SequenceAnnotation {
    /// Length of the annotated span.
    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    /// Whether the span is empty.
    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }
}

/// A composite component — an ordered assembly of parts.
#[derive(Debug, Clone, PartialEq)]
pub struct Component {
    /// Component identifier.
    pub id: String,
    /// Ordered constituent parts (5' → 3').
    pub parts: Vec<Part>,
}

impl Component {
    /// An empty component with the given id.
    pub fn new(id: impl Into<String>) -> Self {
        Component {
            id: id.into(),
            parts: Vec::new(),
        }
    }

    /// Append a part to the 3' end, returning `self` for chaining.
    pub fn with_part(mut self, part: Part) -> Self {
        self.parts.push(part);
        self
    }

    /// The concatenated DNA sequence of every constituent part.
    pub fn sequence(&self) -> Result<Seq> {
        let mut bytes: Vec<u8> = Vec::new();
        for p in &self.parts {
            bytes.extend_from_slice(p.sequence.as_bytes());
        }
        Seq::new(SeqKind::Dna, bytes)
            .map_err(|e| SysbioError::invalid_model("sbol", format!("bad assembly: {e}")))
    }

    /// Total assembled length.
    pub fn len(&self) -> usize {
        self.parts.iter().map(|p| p.len()).sum()
    }

    /// Whether the component has no parts.
    pub fn is_empty(&self) -> bool {
        self.parts.is_empty()
    }

    /// Derive a [`SequenceAnnotation`] for every part, with offsets
    /// computed from the parts' lengths — the SBOL `SequenceAnnotation`
    /// list of the composite component.
    pub fn annotations(&self) -> Vec<SequenceAnnotation> {
        let mut out = Vec::with_capacity(self.parts.len());
        let mut offset = 0;
        for p in &self.parts {
            let end = offset + p.len();
            out.push(SequenceAnnotation {
                id: format!("{}_anno", p.id),
                start: offset,
                end,
                role: p.role,
            });
            offset = end;
        }
        out
    }

    /// Indices of the parts that play a given role.
    pub fn parts_with_role(&self, role: PartRole) -> Vec<usize> {
        self.parts
            .iter()
            .enumerate()
            .filter(|(_, p)| p.role == role)
            .map(|(i, _)| i)
            .collect()
    }
}

/// A functional interaction between two named entities in a device —
/// the SBOL `Interaction` (a `Module`-level relationship).
#[derive(Debug, Clone, PartialEq)]
pub struct Module {
    /// Interaction identifier.
    pub id: String,
    /// The acting / source entity (e.g. a CDS product).
    pub source: String,
    /// The affected / target entity (e.g. a promoter).
    pub target: String,
    /// `true` for activation, `false` for repression.
    pub activates: bool,
}

/// A complete genetic device — components plus their functional
/// interactions. The SBOL `ModuleDefinition`.
#[derive(Debug, Clone, PartialEq)]
pub struct Device {
    /// Device identifier.
    pub id: String,
    /// Structural components.
    pub components: Vec<Component>,
    /// Functional interactions between components / products.
    pub modules: Vec<Module>,
}

impl Device {
    /// An empty device with the given id.
    pub fn new(id: impl Into<String>) -> Self {
        Device {
            id: id.into(),
            components: Vec::new(),
            modules: Vec::new(),
        }
    }

    /// Append a component.
    pub fn add_component(&mut self, c: Component) {
        self.components.push(c);
    }

    /// Append a functional-interaction module.
    pub fn add_module(&mut self, m: Module) {
        self.modules.push(m);
    }

    /// Structural validation: every module endpoint must name a known
    /// component or a part within one.
    pub fn validate(&self) -> Result<()> {
        let known: Vec<&str> = self
            .components
            .iter()
            .flat_map(|c| {
                std::iter::once(c.id.as_str()).chain(c.parts.iter().map(|p| p.id.as_str()))
            })
            .collect();
        for m in &self.modules {
            for endpoint in [&m.source, &m.target] {
                if !known.iter().any(|k| k == endpoint) {
                    return Err(SysbioError::invalid_model(
                        "sbol",
                        format!("module `{}` references unknown entity `{endpoint}`", m.id),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Total assembled length of every component.
    pub fn total_length(&self) -> usize {
        self.components.iter().map(|c| c.len()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn promoter() -> Part {
        Part::new("pTet", PartRole::Promoter, "TCCCTATCAGTGATAGAGA").unwrap()
    }
    fn rbs() -> Part {
        Part::new("B0034", PartRole::Rbs, "AAAGAGGAGAAA").unwrap()
    }
    fn cds() -> Part {
        Part::new("gfp", PartRole::Cds, "ATGGTGAGCAAGGGCGAGGAGTAA").unwrap()
    }

    #[test]
    fn part_carries_role_and_sequence() {
        let p = promoter();
        assert_eq!(p.role, PartRole::Promoter);
        assert_eq!(p.len(), 19);
        assert_eq!(p.role.so_tag(), "SO:promoter");
    }

    #[test]
    fn component_concatenates_part_sequences() {
        let comp = Component::new("tu1")
            .with_part(promoter())
            .with_part(rbs())
            .with_part(cds());
        let seq = comp.sequence().unwrap();
        assert_eq!(seq.len(), 19 + 12 + 24);
        assert_eq!(comp.len(), seq.len());
    }

    #[test]
    fn annotations_have_contiguous_offsets() {
        let comp = Component::new("tu1")
            .with_part(promoter())
            .with_part(rbs())
            .with_part(cds());
        let annos = comp.annotations();
        assert_eq!(annos.len(), 3);
        assert_eq!(annos[0].start, 0);
        assert_eq!(annos[0].end, 19);
        assert_eq!(annos[1].start, 19);
        assert_eq!(annos[2].end, 19 + 12 + 24);
        // No gaps between consecutive annotations.
        for w in annos.windows(2) {
            assert_eq!(w[0].end, w[1].start);
        }
    }

    #[test]
    fn parts_with_role_filters() {
        let comp = Component::new("tu1")
            .with_part(promoter())
            .with_part(rbs())
            .with_part(cds());
        assert_eq!(comp.parts_with_role(PartRole::Cds), vec![2]);
        assert_eq!(
            comp.parts_with_role(PartRole::Terminator),
            Vec::<usize>::new()
        );
    }

    #[test]
    fn device_validates_module_endpoints() {
        let comp = Component::new("tu1").with_part(promoter()).with_part(cds());
        let mut dev = Device::new("d1");
        dev.add_component(comp);
        dev.add_module(Module {
            id: "i1".into(),
            source: "gfp".into(),
            target: "pTet".into(),
            activates: false,
        });
        assert!(dev.validate().is_ok());
    }

    #[test]
    fn device_rejects_dangling_module() {
        let mut dev = Device::new("d1");
        dev.add_component(Component::new("tu1").with_part(promoter()));
        dev.add_module(Module {
            id: "bad".into(),
            source: "ghost".into(),
            target: "pTet".into(),
            activates: true,
        });
        assert!(dev.validate().is_err());
    }

    #[test]
    fn bad_sequence_is_rejected() {
        // 'Z' is not a valid nucleotide.
        assert!(Part::new("x", PartRole::Cds, "ATGZZZ").is_err());
    }
}
