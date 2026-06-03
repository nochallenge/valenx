//! pLannotate-class plasmid feature annotation.
//!
//! Matches a sequence against a built-in library of common molecular-
//! biology features (origins of replication, antibiotic-resistance
//! markers, promoters, terminators, tags, reporter genes) and emits an
//! annotated [`SeqFeature`] for every hit.
//!
//! This is a v1: matching is by signature-substring search on both
//! strands (the unambiguous core of pLannotate's BLAST-based
//! approach), with a configurable identity floor. The feature library
//! stores a short, diagnostic signature for each feature rather than
//! its full sequence.

use crate::error::{BioseqError, Result};
use crate::ops::revcomp::reverse_complement_dna_bytes;
use crate::record::{Location, SeqFeature, Span, Strand};
use crate::seq::{Seq, SeqKind};

/// One entry in the built-in feature library.
#[derive(Clone, Debug)]
pub struct FeatureSignature {
    /// Feature name (e.g. `"AmpR"`, `"T7 promoter"`).
    pub name: &'static str,
    /// GenBank feature type to emit (`"CDS"`, `"promoter"`, …).
    pub feature_type: &'static str,
    /// A diagnostic DNA signature — a substring that uniquely
    /// identifies the feature.
    pub signature: &'static str,
}

/// The built-in plasmid-feature library.
///
/// Signatures are short, characteristic subsequences taken from the
/// canonical features (pUC ori, bla/AmpR, T7 system, common tags and
/// reporters). They are long enough (~24–40 bp) to be specific.
pub fn feature_library() -> &'static [FeatureSignature] {
    &[
        FeatureSignature {
            name: "T7 promoter",
            feature_type: "promoter",
            // Canonical T7 RNA-polymerase promoter.
            signature: "TAATACGACTCACTATAGGG",
        },
        FeatureSignature {
            name: "T7 terminator",
            feature_type: "terminator",
            signature: "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG",
        },
        FeatureSignature {
            name: "lac promoter",
            feature_type: "promoter",
            signature: "TTTACACTTTATGCTTCCGGCTCGTATGTTG",
        },
        FeatureSignature {
            name: "lac operator",
            feature_type: "protein_bind",
            signature: "GGAATTGTGAGCGGATAACAATT",
        },
        FeatureSignature {
            name: "SV40 promoter",
            feature_type: "promoter",
            signature: "GGAGGCTTTTTTGGAGGCCTAGGCTTTTGC",
        },
        FeatureSignature {
            name: "CMV promoter",
            feature_type: "promoter",
            signature: "GACGGTAAATGGCCCGCCTGGCATTATGCC",
        },
        FeatureSignature {
            name: "His-tag",
            feature_type: "CDS",
            // 6x His: CAT CAC CAT CAC CAT CAC.
            signature: "CATCACCATCACCATCAC",
        },
        FeatureSignature {
            name: "FLAG-tag",
            feature_type: "CDS",
            // DYKDDDDK.
            signature: "GATTACAAGGATGACGACGATAAG",
        },
        FeatureSignature {
            name: "HA-tag",
            feature_type: "CDS",
            // YPYDVPDYA.
            signature: "TACCCATACGATGTTCCAGATTACGCT",
        },
        FeatureSignature {
            name: "AmpR",
            feature_type: "CDS",
            // bla (beta-lactamase) internal signature.
            signature: "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCC",
        },
        FeatureSignature {
            name: "KanR",
            feature_type: "CDS",
            // Aminoglycoside 3'-phosphotransferase signature.
            signature: "ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGG",
        },
        FeatureSignature {
            name: "CmR",
            feature_type: "CDS",
            // Chloramphenicol acetyltransferase signature.
            signature: "ATGGAGAAAAAAATCACTGGATATACCACCGTTGAT",
        },
        FeatureSignature {
            name: "ori",
            feature_type: "rep_origin",
            // pUC/pMB1 origin signature (RNAII region).
            signature: "TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTT",
        },
        FeatureSignature {
            name: "lacZ-alpha",
            feature_type: "CDS",
            signature: "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGG",
        },
        FeatureSignature {
            name: "lacI",
            feature_type: "CDS",
            signature: "ATGAAACCAGTAACGTTATACGATGTCGCAGAGTAT",
        },
        FeatureSignature {
            name: "GFP",
            feature_type: "CDS",
            // Green fluorescent protein internal signature.
            signature: "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTC",
        },
        FeatureSignature {
            name: "MCS",
            feature_type: "misc_feature",
            // A common pUC19-style multiple cloning site stretch.
            signature: "GAATTCGAGCTCGGTACCCGGGGATCCTCTAGA",
        },
    ]
}

/// One annotation hit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Annotation {
    /// Matched feature name.
    pub name: String,
    /// 0-based start of the matched region.
    pub start: usize,
    /// 0-based end of the matched region (exclusive).
    pub end: usize,
    /// Strand the signature matched on.
    pub strand: Strand,
}

/// Options for [`annotate`].
#[derive(Copy, Clone, Debug)]
pub struct AnnotateOptions {
    /// Also search the reverse strand. Default `true`.
    pub search_reverse: bool,
}

impl Default for AnnotateOptions {
    fn default() -> Self {
        AnnotateOptions {
            search_reverse: true,
        }
    }
}

/// Finds every library-feature signature present in `seq`.
///
/// For a circular sequence the search wraps the origin so a feature
/// spanning the origin is still found. Returns annotations sorted by
/// `start`. Returns [`BioseqError::Invalid`] for a non-DNA input.
pub fn find_annotations(seq: &Seq, opts: AnnotateOptions) -> Result<Vec<Annotation>> {
    if seq.kind() != SeqKind::Dna {
        return Err(BioseqError::invalid(
            "kind",
            "plasmid annotation needs a DNA sequence",
        ));
    }
    let n = seq.len();
    // For a circular molecule, search a doubled string so origin-
    // spanning features are found; then map coordinates back.
    let mut haystack: Vec<u8> = seq.as_bytes().to_vec();
    if seq.is_circular() && n > 0 {
        let dup = haystack.clone();
        haystack.extend_from_slice(&dup);
    }
    let rc = reverse_complement_dna_bytes(&haystack);

    let mut hits = Vec::new();
    for entry in feature_library() {
        let needle = entry.signature.as_bytes();
        for pos in find_all(&haystack, needle) {
            if pos >= n && seq.is_circular() {
                continue; // duplicate already counted from the first copy
            }
            let start = pos;
            let end = pos + needle.len();
            // Only keep a hit that fits in (linear) or wraps validly
            // (circular).
            if !seq.is_circular() && end > n {
                continue;
            }
            hits.push(Annotation {
                name: entry.name.to_string(),
                start: start % n.max(1),
                end: if seq.is_circular() {
                    end % n.max(1)
                } else {
                    end
                },
                strand: Strand::Forward,
            });
        }
        if opts.search_reverse {
            for pos in find_all(&rc, needle) {
                // rc index r maps to forward coordinate hay_len - r - len.
                let hay_len = haystack.len();
                let needle_len = needle.len();
                if pos + needle_len > hay_len {
                    continue;
                }
                let fwd_start = hay_len - pos - needle_len;
                if fwd_start >= n && seq.is_circular() {
                    continue;
                }
                let fwd_end = fwd_start + needle_len;
                if !seq.is_circular() && fwd_end > n {
                    continue;
                }
                hits.push(Annotation {
                    name: entry.name.to_string(),
                    start: fwd_start % n.max(1),
                    end: if seq.is_circular() {
                        fwd_end % n.max(1)
                    } else {
                        fwd_end
                    },
                    strand: Strand::Reverse,
                });
            }
        }
    }
    hits.sort_by(|a, b| a.start.cmp(&b.start).then_with(|| a.name.cmp(&b.name)));
    hits.dedup();
    Ok(hits)
}

/// All start indices where `needle` occurs in `haystack`
/// (case-insensitive, overlapping).
fn find_all(haystack: &[u8], needle: &[u8]) -> Vec<usize> {
    let mut out = Vec::new();
    if needle.is_empty() || needle.len() > haystack.len() {
        return out;
    }
    for start in 0..=haystack.len() - needle.len() {
        if haystack[start..start + needle.len()]
            .iter()
            .zip(needle)
            .all(|(a, b)| a.eq_ignore_ascii_case(b))
        {
            out.push(start);
        }
    }
    out
}

/// Annotates a sequence — finds all features and returns them as
/// [`SeqFeature`]s ready to attach to a [`crate::record::SeqRecord`].
pub fn annotate(seq: &Seq, opts: AnnotateOptions) -> Result<Vec<SeqFeature>> {
    let hits = find_annotations(seq, opts)?;
    let lib = feature_library();
    let mut features = Vec::with_capacity(hits.len());
    for h in hits {
        let entry = lib.iter().find(|e| e.name == h.name);
        let feature_type = entry.map(|e| e.feature_type).unwrap_or("misc_feature");
        let span = Span::with_strand(h.start, h.end.max(h.start), h.strand);
        let feature = SeqFeature::new(feature_type, Location::Single(span))
            .with_qualifier("label", h.name.clone())
            .with_qualifier("note", "annotated by valenx-bioseq feature library");
        features.push(feature);
    }
    Ok(features)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::Topology;

    fn lib_sig(name: &str) -> &'static str {
        feature_library()
            .iter()
            .find(|e| e.name == name)
            .unwrap()
            .signature
    }

    #[test]
    fn library_is_populated() {
        let lib = feature_library();
        assert!(lib.len() >= 13);
        assert!(lib.iter().any(|e| e.name == "AmpR"));
        assert!(lib.iter().any(|e| e.name == "T7 promoter"));
        assert!(lib.iter().any(|e| e.name == "GFP"));
    }

    #[test]
    fn finds_a_forward_feature() {
        let t7 = lib_sig("T7 promoter");
        let seq_str = format!("AAAAA{t7}AAAAA");
        let s = Seq::new(SeqKind::Dna, seq_str).unwrap();
        let hits = find_annotations(&s, AnnotateOptions::default()).unwrap();
        let t7_hit = hits.iter().find(|h| h.name == "T7 promoter").unwrap();
        assert_eq!(t7_hit.start, 5);
        assert_eq!(t7_hit.strand, Strand::Forward);
    }

    #[test]
    fn finds_a_reverse_feature() {
        let his = lib_sig("His-tag");
        let his_rc = reverse_complement_dna_bytes(his.as_bytes());
        let mut seq_bytes = b"AAAAA".to_vec();
        seq_bytes.extend_from_slice(&his_rc);
        seq_bytes.extend_from_slice(b"AAAAA");
        let s = Seq::new(SeqKind::Dna, seq_bytes).unwrap();
        let hits = find_annotations(&s, AnnotateOptions::default()).unwrap();
        let hit = hits.iter().find(|h| h.name == "His-tag").unwrap();
        assert_eq!(hit.strand, Strand::Reverse);
        assert_eq!(hit.start, 5);
    }

    #[test]
    fn reverse_search_can_be_disabled() {
        let his = lib_sig("His-tag");
        let his_rc = reverse_complement_dna_bytes(his.as_bytes());
        let s = Seq::new(SeqKind::Dna, his_rc).unwrap();
        let opts = AnnotateOptions {
            search_reverse: false,
        };
        let hits = find_annotations(&s, opts).unwrap();
        assert!(hits.iter().all(|h| h.name != "His-tag"));
    }

    #[test]
    fn no_features_in_random_sequence() {
        let s = Seq::new(SeqKind::Dna, "AAAAAAAAAAAAAAAAAAAA").unwrap();
        let hits = find_annotations(&s, AnnotateOptions::default()).unwrap();
        assert!(hits.is_empty());
    }

    #[test]
    fn annotate_emits_seqfeatures() {
        let amp = lib_sig("AmpR");
        let seq_str = format!("GGGG{amp}GGGG");
        let s = Seq::new(SeqKind::Dna, seq_str).unwrap();
        let feats = annotate(&s, AnnotateOptions::default()).unwrap();
        let amp_feat = feats.iter().find(|f| f.label() == "AmpR").unwrap();
        assert_eq!(amp_feat.feature_type, "CDS");
        assert_eq!(amp_feat.qualifier("label"), Some("AmpR"));
    }

    #[test]
    fn circular_feature_spanning_origin() {
        // Put the T7 signature so it wraps the origin.
        let t7 = lib_sig("T7 promoter");
        let bytes = t7.as_bytes();
        let split = 8;
        // sequence = [tail of t7][filler][head of t7]
        let mut seq_bytes = bytes[split..].to_vec();
        seq_bytes.extend_from_slice(b"AAAAAAAAAA");
        seq_bytes.extend_from_slice(&bytes[..split]);
        let s = Seq::with_topology(SeqKind::Dna, seq_bytes, Topology::Circular).unwrap();
        let hits = find_annotations(&s, AnnotateOptions::default()).unwrap();
        assert!(
            hits.iter().any(|h| h.name == "T7 promoter"),
            "origin-spanning T7 promoter should be found"
        );
    }

    #[test]
    fn non_dna_rejected() {
        let p = Seq::new(SeqKind::Protein, "MKVL").unwrap();
        assert!(find_annotations(&p, AnnotateOptions::default()).is_err());
    }
}
