//! Feature 5 — multiplex editing and gRNA-array design.
//!
//! Editing several loci at once needs several guides delivered
//! together. Two delivery formats dominate:
//!
//! - a **Cas9 multi-cassette array** — each guide gets its own Pol-III
//!   (U6 / H1) promoter and tracrRNA scaffold;
//! - a **Cas12a crRNA array** — a single transcript carrying every
//!   guide separated by the Cas12a direct-repeat, which Cas12a
//!   processes itself into individual crRNAs.
//!
//! This module batch-designs guides across many targets
//! ([`crate::crispr::guide_design`]), picks the best guide per target,
//! checks the chosen guides for **cross-target collisions** (two
//! guides whose protospacers are near-identical can mis-prime), and
//! assembles a gRNA-array sequence.
//!
//! ## v1 scope
//!
//! Array assembly emits a representative direct-repeat / scaffold
//! layout, not a vendor-specific cloning construct. The collision
//! check is a Hamming-distance screen between chosen protospacers; it
//! does not run a full off-target scan of every guide against every
//! other target's flanks.

use crate::crispr::guide_design::{design_guides, GuideCandidate, GuideDesignRequest};
use crate::crispr::nuclease::{NucleaseClass, NucleaseId};
use crate::error::{GeneditingError, Result};
use serde::{Deserialize, Serialize};

/// One target locus in a multiplex design.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MultiplexTarget {
    /// A human-readable label (gene name, locus id).
    pub label: String,
    /// The target region (forward-strand ACGT window).
    pub region: Vec<u8>,
}

/// A request for a multiplex / gRNA-array design.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MultiplexRequest {
    /// The loci to edit.
    pub targets: Vec<MultiplexTarget>,
    /// Which nuclease — drives both guide design and the array format
    /// (Cas12a ⇒ crRNA array, Cas9 ⇒ multi-cassette).
    pub nuclease: NucleaseId,
    /// Minimum Hamming distance required between any two chosen
    /// protospacers; closer pairs are flagged as a collision risk.
    pub min_guide_distance: usize,
}

impl MultiplexRequest {
    /// A request with a sensible default collision threshold (6).
    pub fn new(targets: Vec<MultiplexTarget>, nuclease: NucleaseId) -> Self {
        MultiplexRequest {
            targets,
            nuclease,
            min_guide_distance: 6,
        }
    }
}

/// The chosen guide for one target in a multiplex design.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MultiplexGuide {
    /// The target label this guide edits.
    pub label: String,
    /// The chosen guide candidate.
    pub guide: GuideCandidate,
}

/// A collision between two chosen guides (a mis-priming risk).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GuideCollision {
    /// Label of the first guide's target.
    pub label_a: String,
    /// Label of the second guide's target.
    pub label_b: String,
    /// Hamming distance between the two protospacers.
    pub distance: usize,
}

/// The result of a multiplex / gRNA-array design.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MultiplexDesign {
    /// One chosen guide per successfully designed target.
    pub guides: Vec<MultiplexGuide>,
    /// Targets for which no guide could be designed, with the reason.
    pub failed_targets: Vec<(String, String)>,
    /// Cross-target protospacer collisions (mis-priming risk).
    pub collisions: Vec<GuideCollision>,
    /// The assembled gRNA-array transcript / construct sequence.
    pub array: Vec<u8>,
    /// A description of the array format used.
    pub array_format: String,
}

impl MultiplexDesign {
    /// `true` when no cross-target collision was detected.
    pub fn is_collision_free(&self) -> bool {
        self.collisions.is_empty()
    }
}

/// Hamming distance between two equal-length byte slices; for unequal
/// lengths every position past the shorter one counts as a mismatch.
fn hamming(a: &[u8], b: &[u8]) -> usize {
    let common = a.len().min(b.len());
    let mut d = a.len().abs_diff(b.len());
    for i in 0..common {
        if !a[i].eq_ignore_ascii_case(&b[i]) {
            d += 1;
        }
    }
    d
}

/// The Cas12a direct repeat (the wild-type AsCas12a 19 nt repeat) used
/// to separate guides in a self-processing crRNA array.
const CAS12A_DIRECT_REPEAT: &[u8] = b"AATTTCTACTAAGTGTAGAT";

/// A representative U6-promoter + tracrRNA-scaffold spacer used between
/// guides in a Cas9 multi-cassette array (shortened sentinel — the
/// real cassette is far longer; this marks the boundary).
const CAS9_CASSETTE_LINKER: &[u8] = b"TTTTTTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC";

/// Designs a multiplex edit and assembles a gRNA array.
///
/// Designs guides for every target, keeps the single best guide per
/// target, screens the chosen guides for cross-target collisions, and
/// assembles an array transcript appropriate to the nuclease (Cas12a:
/// direct-repeat-separated crRNA array; Cas9: scaffold-linked
/// multi-cassette).
///
/// # Errors
/// - [`GeneditingError::Invalid`] if the target list is empty.
/// - [`GeneditingError::NoValidDesign`] if *no* target yields a guide
///   (per-target failures are otherwise collected, not fatal).
pub fn design_multiplex(req: &MultiplexRequest) -> Result<MultiplexDesign> {
    if req.targets.is_empty() {
        return Err(GeneditingError::invalid(
            "targets",
            "multiplex design needs at least one target",
        ));
    }
    let nuc = crate::crispr::nuclease::nuclease(req.nuclease);

    let mut guides: Vec<MultiplexGuide> = Vec::new();
    let mut failed: Vec<(String, String)> = Vec::new();

    for t in &req.targets {
        let g_req = GuideDesignRequest::new(t.region.clone(), req.nuclease);
        match design_guides(&g_req) {
            Ok(cands) => {
                // design_guides already sorts best-first.
                guides.push(MultiplexGuide {
                    label: t.label.clone(),
                    guide: cands.into_iter().next().expect("non-empty on Ok"),
                });
            }
            Err(e) => failed.push((t.label.clone(), e.to_string())),
        }
    }
    if guides.is_empty() {
        return Err(GeneditingError::no_valid_design(
            "multiplex",
            "no target in the set yields a usable guide",
        ));
    }

    // Cross-target collision screen.
    let mut collisions = Vec::new();
    for i in 0..guides.len() {
        for j in (i + 1)..guides.len() {
            let d = hamming(
                guides[i].guide.protospacer.as_bytes(),
                guides[j].guide.protospacer.as_bytes(),
            );
            if d < req.min_guide_distance {
                collisions.push(GuideCollision {
                    label_a: guides[i].label.clone(),
                    label_b: guides[j].label.clone(),
                    distance: d,
                });
            }
        }
    }

    // Array assembly.
    let (array, array_format) = assemble_array(&guides, nuc.class);

    Ok(MultiplexDesign {
        guides,
        failed_targets: failed,
        collisions,
        array,
        array_format,
    })
}

/// Assembles the gRNA-array sequence for the chosen guides.
fn assemble_array(guides: &[MultiplexGuide], class: NucleaseClass) -> (Vec<u8>, String) {
    let mut array = Vec::new();
    match class {
        NucleaseClass::Cas12 => {
            // crRNA array: [DR][spacer][DR][spacer]...[DR]
            for g in guides {
                array.extend_from_slice(CAS12A_DIRECT_REPEAT);
                array.extend_from_slice(g.guide.protospacer.as_bytes());
            }
            array.extend_from_slice(CAS12A_DIRECT_REPEAT);
            (
                array,
                format!(
                    "Cas12a self-processing crRNA array — {} spacers separated by the \
                     AsCas12a direct repeat",
                    guides.len()
                ),
            )
        }
        // Cas9 / Cas13: a multi-cassette layout (one scaffold per guide).
        _ => {
            for (i, g) in guides.iter().enumerate() {
                if i > 0 {
                    array.extend_from_slice(CAS9_CASSETTE_LINKER);
                }
                array.extend_from_slice(g.guide.protospacer.as_bytes());
                array.extend_from_slice(CAS9_CASSETTE_LINKER);
            }
            (
                array,
                format!(
                    "Cas9 multi-cassette array — {} guides, each with its own \
                     Pol-III promoter + tracrRNA scaffold (linker shown as a sentinel)",
                    guides.len()
                ),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn targets() -> Vec<MultiplexTarget> {
        vec![
            MultiplexTarget {
                label: "GENE_A".to_string(),
                region: b"ACGTACGTACGTACGTACGTAGGTTTT".to_vec(),
            },
            MultiplexTarget {
                label: "GENE_B".to_string(),
                region: b"TTGGCATCAGTACCGATGCACTGCGGTT".to_vec(),
            },
        ]
    }

    #[test]
    fn designs_guides_for_each_target() {
        let req = MultiplexRequest::new(targets(), NucleaseId::SpCas9);
        let d = design_multiplex(&req).unwrap();
        assert_eq!(d.guides.len(), 2);
        assert!(d.guides.iter().any(|g| g.label == "GENE_A"));
        assert!(d.guides.iter().any(|g| g.label == "GENE_B"));
    }

    #[test]
    fn empty_target_list_is_rejected() {
        let req = MultiplexRequest::new(Vec::new(), NucleaseId::SpCas9);
        assert_eq!(design_multiplex(&req).unwrap_err().code(), "genediting.invalid");
    }

    #[test]
    fn cas9_array_uses_multi_cassette() {
        let req = MultiplexRequest::new(targets(), NucleaseId::SpCas9);
        let d = design_multiplex(&req).unwrap();
        assert!(d.array_format.contains("multi-cassette"));
        assert!(!d.array.is_empty());
    }

    #[test]
    fn cas12a_array_uses_direct_repeats() {
        let t = vec![MultiplexTarget {
            label: "L1".to_string(),
            region: b"TTTAACGTACGTACGTACGTACGTACG".to_vec(),
        }];
        let req = MultiplexRequest::new(t, NucleaseId::Cas12a);
        let d = design_multiplex(&req).unwrap();
        assert!(d.array_format.contains("crRNA array"));
        // The array begins and ends with the direct repeat.
        assert!(d.array.starts_with(CAS12A_DIRECT_REPEAT));
        assert!(d.array.ends_with(CAS12A_DIRECT_REPEAT));
    }

    #[test]
    fn identical_targets_collide() {
        // Two copies of the same locus → identical chosen guides →
        // Hamming distance 0 → a collision.
        let same = b"ACGTACGTACGTACGTACGTAGGTTTT".to_vec();
        let t = vec![
            MultiplexTarget {
                label: "X".to_string(),
                region: same.clone(),
            },
            MultiplexTarget {
                label: "Y".to_string(),
                region: same,
            },
        ];
        let req = MultiplexRequest::new(t, NucleaseId::SpCas9);
        let d = design_multiplex(&req).unwrap();
        assert!(!d.is_collision_free());
        assert_eq!(d.collisions[0].distance, 0);
    }

    #[test]
    fn distinct_guides_are_collision_free() {
        let req = MultiplexRequest::new(targets(), NucleaseId::SpCas9);
        let d = design_multiplex(&req).unwrap();
        assert!(d.is_collision_free());
    }

    #[test]
    fn failed_target_is_collected_not_fatal() {
        let mut t = targets();
        // A target with no PAM site at all.
        t.push(MultiplexTarget {
            label: "NOPAM".to_string(),
            region: b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec(),
        });
        let req = MultiplexRequest::new(t, NucleaseId::SpCas9);
        let d = design_multiplex(&req).unwrap();
        assert!(d.failed_targets.iter().any(|(l, _)| l == "NOPAM"));
        assert_eq!(d.guides.len(), 2);
    }

    #[test]
    fn hamming_helper() {
        assert_eq!(hamming(b"ACGT", b"ACGT"), 0);
        assert_eq!(hamming(b"ACGT", b"ACGA"), 1);
        assert_eq!(hamming(b"ACGT", b"AC"), 2);
    }
}
