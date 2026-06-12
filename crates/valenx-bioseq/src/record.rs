//! [`SeqRecord`] / [`SeqFeature`] / [`Location`] — the annotated-sequence
//! data model.
//!
//! This is the Biopython `SeqRecord` / `SeqFeature` analogue. A
//! [`SeqRecord`] bundles an identifier, a name, a free-text
//! description, the [`Seq`] itself, and a list of [`SeqFeature`]
//! annotations. Each feature has a type string, a [`Location`]
//! (possibly a `join(...)` / `complement(...)` compound location), and
//! a qualifiers map.

use crate::seq::Seq;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Which strand a feature lies on.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize, Default)]
pub enum Strand {
    /// The given (top, 5′→3′) strand.
    #[default]
    Forward,
    /// The reverse-complement (bottom) strand.
    Reverse,
    /// Strand unknown / not applicable.
    Unknown,
}

impl Strand {
    /// `-1` / `+1` / `0` — the GenBank-style numeric strand.
    pub fn sign(self) -> i8 {
        match self {
            Strand::Forward => 1,
            Strand::Reverse => -1,
            Strand::Unknown => 0,
        }
    }

    /// The opposite strand (`Unknown` stays `Unknown`).
    pub fn flip(self) -> Strand {
        match self {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
            Strand::Unknown => Strand::Unknown,
        }
    }
}

/// A contiguous half-open `[start, end)` interval on one strand.
///
/// Coordinates are 0-based — note GenBank flat files are 1-based and
/// inclusive; the GenBank reader/writer converts.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct Span {
    /// 0-based start (inclusive).
    pub start: usize,
    /// 0-based end (exclusive).
    pub end: usize,
    /// Strand the span lies on.
    pub strand: Strand,
}

impl Span {
    /// Builds a forward-strand span.
    pub fn new(start: usize, end: usize) -> Self {
        Span {
            start,
            end,
            strand: Strand::Forward,
        }
    }

    /// Builds a span on an explicit strand.
    pub fn with_strand(start: usize, end: usize, strand: Strand) -> Self {
        Span { start, end, strand }
    }

    /// Length of the span in residues.
    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    /// `true` if the span covers zero residues.
    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }
}

/// A feature location — either a single [`Span`], a compound
/// `join(...)` / `order(...)` of several spans (for multi-exon genes
/// or grouped sites), or a `^`-between-bases position.
///
/// The distinction between [`Location::Join`] and [`Location::Order`]
/// is semantic: `join(...)` asserts that the segments form one
/// contiguous biological entity (e.g. the exons of a spliced mRNA),
/// while `order(...)` only asserts that the segments are listed in a
/// particular order without claiming they are joined (e.g. a primer
/// pair, scattered binding sites). Both are extracted identically.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Location {
    /// A single contiguous interval.
    Single(Span),
    /// `join(...)` — an ordered list of spans making up one feature
    /// (e.g. the exons of a spliced CDS).
    Join(Vec<Span>),
    /// `order(...)` — an ordered list of spans with no joining
    /// claim. Same span data as `Join`; different semantics.
    Order(Vec<Span>),
    /// `n^n+1` — the bond between two adjacent bases.
    /// `Between(n)` refers to the phosphodiester bond between
    /// 0-based positions `n-1` and `n` (the same convention as the
    /// crate's half-open intervals: a `Between(n)` lies *before* base
    /// `n`).
    ///
    /// Used in INSDC records for cut sites, junction points, and other
    /// inter-base annotations (e.g. `/note="cleavage between residues
    /// 100^101"`).
    Between {
        /// The 0-based position the bond lies before.
        position: usize,
        /// Which strand the bond is annotated on.
        strand: Strand,
    },
}

impl Location {
    /// Builds a forward single-span location.
    pub fn single(start: usize, end: usize) -> Self {
        Location::Single(Span::new(start, end))
    }

    /// Builds a between-bases location at the bond before 0-based
    /// position `position`.
    pub fn between(position: usize) -> Self {
        Location::Between {
            position,
            strand: Strand::Forward,
        }
    }

    /// All spans, in order, regardless of variant. A
    /// [`Location::Between`] yields a zero-length span at the bond
    /// position (so callers that just walk spans still see something
    /// sensible).
    pub fn spans(&self) -> Vec<Span> {
        match self {
            Location::Single(s) => vec![*s],
            Location::Join(v) | Location::Order(v) => v.clone(),
            Location::Between { position, strand } => {
                vec![Span::with_strand(*position, *position, *strand)]
            }
        }
    }

    /// Total residue length summed across all spans. Zero for a
    /// [`Location::Between`].
    pub fn total_len(&self) -> usize {
        self.spans().iter().map(Span::len).sum()
    }

    /// The minimum start coordinate across all spans.
    pub fn min_start(&self) -> usize {
        self.spans().iter().map(|s| s.start).min().unwrap_or(0)
    }

    /// The maximum end coordinate across all spans.
    pub fn max_end(&self) -> usize {
        self.spans().iter().map(|s| s.end).max().unwrap_or(0)
    }

    /// `true` if any span is on the reverse strand.
    pub fn is_reverse(&self) -> bool {
        match self {
            Location::Between { strand, .. } => *strand == Strand::Reverse,
            _ => self.spans().iter().any(|s| s.strand == Strand::Reverse),
        }
    }

    /// Extracts the feature's residues from `parent`, honoring strand
    /// and joining the spans. Reverse-strand spans are
    /// reverse-complemented (for nucleotide parents). Out-of-range
    /// spans contribute nothing. A [`Location::Between`] always
    /// extracts an empty sequence (a bond has no residues).
    pub fn extract(&self, parent: &Seq) -> Seq {
        if matches!(self, Location::Between { .. }) {
            return Seq::new_unchecked(parent.kind(), Vec::new(), crate::seq::Topology::Linear);
        }
        use crate::ops::revcomp::reverse_complement;
        let mut out: Vec<u8> = Vec::new();
        let mut spans = self.spans();
        // For a reverse-strand join, GenBank lists exons in reverse
        // genomic order; extraction concatenates them as listed.
        let any_rev = self.is_reverse();
        if any_rev {
            // Walk spans; each reverse span contributes its rev-comp.
            for s in &spans {
                let piece = clamp_slice(parent, s.start, s.end);
                if s.strand == Strand::Reverse {
                    if let Ok(rc) = reverse_complement(&piece) {
                        out.extend_from_slice(rc.as_bytes());
                    } else {
                        out.extend_from_slice(piece.as_bytes());
                    }
                } else {
                    out.extend_from_slice(piece.as_bytes());
                }
            }
        } else {
            spans.sort_by_key(|s| s.start);
            for s in &spans {
                let piece = clamp_slice(parent, s.start, s.end);
                out.extend_from_slice(piece.as_bytes());
            }
        }
        Seq::new_unchecked(parent.kind(), out, crate::seq::Topology::Linear)
    }
}

/// Slices `parent[start..end]` clamping to bounds; never panics.
fn clamp_slice(parent: &Seq, start: usize, end: usize) -> Seq {
    let n = parent.len();
    let s = start.min(n);
    let e = end.min(n);
    parent.slice(s, e.max(s)).unwrap_or_else(|_| {
        Seq::new_unchecked(parent.kind(), Vec::new(), crate::seq::Topology::Linear)
    })
}

/// An annotation on a [`SeqRecord`] — the Biopython `SeqFeature`
/// analogue.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SeqFeature {
    /// Feature type — `"CDS"`, `"gene"`, `"promoter"`, `"misc_feature"`…
    pub feature_type: String,
    /// Where the feature lies on the parent sequence.
    pub location: Location,
    /// Free-form qualifiers (`/gene="lacZ"`, `/product="..."` …).
    /// Ordered for deterministic output.
    pub qualifiers: BTreeMap<String, String>,
}

impl SeqFeature {
    /// Builds a feature with an empty qualifier map.
    pub fn new(feature_type: impl Into<String>, location: Location) -> Self {
        SeqFeature {
            feature_type: feature_type.into(),
            location,
            qualifiers: BTreeMap::new(),
        }
    }

    /// Builder-style qualifier insertion.
    pub fn with_qualifier(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.qualifiers.insert(key.into(), value.into());
        self
    }

    /// Looks up a qualifier value.
    pub fn qualifier(&self, key: &str) -> Option<&str> {
        self.qualifiers.get(key).map(String::as_str)
    }

    /// A human-friendly label — the `label`, `gene`, or `product`
    /// qualifier, falling back to the feature type.
    pub fn label(&self) -> &str {
        self.qualifier("label")
            .or_else(|| self.qualifier("gene"))
            .or_else(|| self.qualifier("product"))
            .unwrap_or(&self.feature_type)
    }
}

/// A literature reference attached to a [`SeqRecord`].
///
/// Mirrors the GenBank `REFERENCE` block and the EMBL `RN`/`RC`/`RP`/
/// `RA`/`RT`/`RL`/`RX` line group. Free-text fields are preserved
/// verbatim; only the structural layout is parsed.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Reference {
    /// The reference number (1-based), as listed by `REFERENCE` / `RN`.
    pub number: usize,
    /// The base range the reference covers, e.g. `"1..1859"`. Stored
    /// as the raw flat-file text — many citations specify a list of
    /// ranges and we deliberately do not reformat them.
    pub bases: String,
    /// `AUTHORS` / `RA` — comma-separated author list, verbatim.
    pub authors: String,
    /// `CONSRTM` — consortium name, if present.
    pub consortium: String,
    /// `TITLE` / `RT` — article title.
    pub title: String,
    /// `JOURNAL` / `RL` — journal citation (volume, page, year).
    pub journal: String,
    /// `PUBMED` / `RX PUBMED;...` — PubMed identifier, if present.
    pub pubmed: String,
    /// `MEDLINE` — MEDLINE identifier (legacy), if present.
    pub medline: String,
    /// `REMARK` / `RC` — free-text remark / comment.
    pub remark: String,
}

/// An identified, annotated sequence — the Biopython `SeqRecord`
/// analogue.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SeqRecord {
    /// Primary identifier (accession / locus name).
    pub id: String,
    /// Short name.
    pub name: String,
    /// Free-text description.
    pub description: String,
    /// The sequence itself.
    pub seq: Seq,
    /// Feature annotations.
    pub features: Vec<SeqFeature>,
    /// Free-form record-level annotations (`organism`, `date`, …).
    pub annotations: BTreeMap<String, String>,
    /// Literature references — `REFERENCE` blocks (GenBank) / `RN/RA/
    /// RT/RL/RX` line groups (EMBL).
    #[serde(default)]
    pub references: Vec<Reference>,
}

impl SeqRecord {
    /// Builds a record with the given id and sequence; name and
    /// description default to the id / empty.
    pub fn new(id: impl Into<String>, seq: Seq) -> Self {
        let id = id.into();
        SeqRecord {
            name: id.clone(),
            description: String::new(),
            id,
            seq,
            features: Vec::new(),
            annotations: BTreeMap::new(),
            references: Vec::new(),
        }
    }

    /// Builder-style description setter.
    pub fn with_description(mut self, d: impl Into<String>) -> Self {
        self.description = d.into();
        self
    }

    /// Appends a feature.
    pub fn add_feature(&mut self, feature: SeqFeature) {
        self.features.push(feature);
    }

    /// All features of a given type (e.g. `"CDS"`).
    pub fn features_of_type<'a>(&'a self, ty: &'a str) -> impl Iterator<Item = &'a SeqFeature> {
        self.features.iter().filter(move |f| f.feature_type == ty)
    }

    /// Number of residues in the record's sequence.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// `true` if the record's sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::SeqKind;

    #[test]
    fn span_and_location_geometry() {
        let loc = Location::Join(vec![Span::new(0, 3), Span::new(10, 16)]);
        assert_eq!(loc.total_len(), 9);
        assert_eq!(loc.min_start(), 0);
        assert_eq!(loc.max_end(), 16);
        assert!(!loc.is_reverse());
    }

    #[test]
    fn location_extract_forward_join() {
        let parent = Seq::new(SeqKind::Dna, "AAACCCGGGTTT").unwrap();
        // join the first 3 and last 3 residues
        let loc = Location::Join(vec![Span::new(0, 3), Span::new(9, 12)]);
        assert_eq!(loc.extract(&parent).as_str(), "AAATTT");
    }

    #[test]
    fn location_extract_reverse() {
        let parent = Seq::new(SeqKind::Dna, "AAACCCGGGTTT").unwrap();
        let loc = Location::Single(Span::with_strand(0, 3, Strand::Reverse));
        // AAA reverse-complemented -> TTT
        assert_eq!(loc.extract(&parent).as_str(), "TTT");
    }

    #[test]
    fn extract_clamps_out_of_range() {
        let parent = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        let loc = Location::single(2, 100);
        assert_eq!(loc.extract(&parent).as_str(), "GT");
    }

    #[test]
    fn feature_qualifiers_and_label() {
        let f = SeqFeature::new("CDS", Location::single(0, 30))
            .with_qualifier("gene", "lacZ")
            .with_qualifier("product", "beta-galactosidase");
        assert_eq!(f.qualifier("gene"), Some("lacZ"));
        assert_eq!(f.label(), "lacZ");
        let g = SeqFeature::new("misc_feature", Location::single(0, 3));
        assert_eq!(g.label(), "misc_feature");
    }

    #[test]
    fn record_roundtrip() {
        let seq = Seq::new(SeqKind::Dna, "ATGAAATAG").unwrap();
        let mut rec = SeqRecord::new("test1", seq).with_description("demo");
        rec.add_feature(SeqFeature::new("CDS", Location::single(0, 9)));
        assert_eq!(rec.len(), 9);
        assert_eq!(rec.features_of_type("CDS").count(), 1);
        assert_eq!(rec.description, "demo");
    }

    #[test]
    fn strand_helpers() {
        assert_eq!(Strand::Forward.sign(), 1);
        assert_eq!(Strand::Reverse.sign(), -1);
        assert_eq!(Strand::Forward.flip(), Strand::Reverse);
    }

    #[test]
    fn order_location_extracts_like_join() {
        let parent = Seq::new(SeqKind::Dna, "AAACCCGGGTTT").unwrap();
        let join = Location::Join(vec![Span::new(0, 3), Span::new(9, 12)]);
        let order = Location::Order(vec![Span::new(0, 3), Span::new(9, 12)]);
        // Same residues; the distinction is semantic, not arithmetic.
        assert_eq!(join.extract(&parent).as_str(), "AAATTT");
        assert_eq!(order.extract(&parent).as_str(), "AAATTT");
        assert_eq!(order.total_len(), 6);
        assert_eq!(order.spans().len(), 2);
    }

    #[test]
    fn between_location_geometry_and_extract() {
        let parent = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        let bond = Location::between(2);
        // A bond has zero residues and extracts an empty sequence.
        assert_eq!(bond.total_len(), 0);
        assert_eq!(bond.extract(&parent).as_str(), "");
        // The single "span" sits at the position.
        let spans = bond.spans();
        assert_eq!(spans.len(), 1);
        assert_eq!(spans[0].start, 2);
        assert_eq!(spans[0].end, 2);
        assert!(spans[0].is_empty());
    }
}
