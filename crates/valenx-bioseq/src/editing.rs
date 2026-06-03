//! Sequence editing — insert, delete, replace, and feature-coordinate
//! updates after an edit.
//!
//! The functions here operate on a [`SeqRecord`] (or a bare [`Seq`])
//! and, crucially, keep the record's [`SeqFeature`] annotations in
//! sync: when bases are inserted or deleted, every downstream feature
//! span is shifted, and features overlapping a deletion are clipped or
//! dropped.

use crate::error::{BioseqError, Result};
use crate::record::{Location, SeqFeature, SeqRecord, Span};
use crate::seq::Seq;

/// Inserts `insert` into `seq` at 0-based position `pos`.
///
/// `pos == seq.len()` appends. Returns [`BioseqError::Invalid`] for a
/// kind mismatch or an out-of-range position.
pub fn insert(seq: &Seq, pos: usize, insert: &Seq) -> Result<Seq> {
    if seq.kind() != insert.kind() {
        return Err(BioseqError::invalid("kind", "insert kind mismatch"));
    }
    if pos > seq.len() {
        return Err(BioseqError::invalid(
            "pos",
            format!("position {pos} exceeds length {}", seq.len()),
        ));
    }
    let mut bytes = seq.as_bytes()[..pos].to_vec();
    bytes.extend_from_slice(insert.as_bytes());
    bytes.extend_from_slice(&seq.as_bytes()[pos..]);
    Ok(Seq::new_unchecked(seq.kind(), bytes, seq.topology()))
}

/// Deletes the half-open range `[start, end)` from `seq`.
///
/// Returns [`BioseqError::Invalid`] for an inverted or out-of-range
/// range.
pub fn delete(seq: &Seq, start: usize, end: usize) -> Result<Seq> {
    if start > end || end > seq.len() {
        return Err(BioseqError::invalid(
            "range",
            format!("bad delete range [{start},{end}) on length {}", seq.len()),
        ));
    }
    let mut bytes = seq.as_bytes()[..start].to_vec();
    bytes.extend_from_slice(&seq.as_bytes()[end..]);
    Ok(Seq::new_unchecked(seq.kind(), bytes, seq.topology()))
}

/// Replaces the half-open range `[start, end)` of `seq` with
/// `replacement` (which may be a different length).
pub fn replace(seq: &Seq, start: usize, end: usize, replacement: &Seq) -> Result<Seq> {
    if seq.kind() != replacement.kind() {
        return Err(BioseqError::invalid("kind", "replacement kind mismatch"));
    }
    if start > end || end > seq.len() {
        return Err(BioseqError::invalid(
            "range",
            format!("bad replace range [{start},{end}) on length {}", seq.len()),
        ));
    }
    let mut bytes = seq.as_bytes()[..start].to_vec();
    bytes.extend_from_slice(replacement.as_bytes());
    bytes.extend_from_slice(&seq.as_bytes()[end..]);
    Ok(Seq::new_unchecked(seq.kind(), bytes, seq.topology()))
}

/// Shifts a single [`Span`] to account for an insertion of `len` bases
/// at `pos`. A span entirely before `pos` is unchanged; a span at or
/// after `pos` shifts by `len`; a span straddling `pos` grows by `len`.
fn shift_span_for_insert(span: Span, pos: usize, len: usize) -> Span {
    let new_start = if span.start >= pos {
        span.start + len
    } else {
        span.start
    };
    let new_end = if span.end > pos {
        span.end + len
    } else {
        span.end
    };
    Span {
        start: new_start,
        end: new_end,
        strand: span.strand,
    }
}

/// Adjusts a single [`Span`] for a deletion of `[del_start, del_end)`.
/// Returns `None` if the span is entirely consumed by the deletion.
fn shift_span_for_delete(span: Span, del_start: usize, del_end: usize) -> Option<Span> {
    let del_len = del_end - del_start;
    // Clip the span to outside the deleted region.
    let clip = |coord: usize| -> usize {
        if coord <= del_start {
            coord
        } else if coord >= del_end {
            coord - del_len
        } else {
            // Inside the deletion -> snaps to the deletion start.
            del_start
        }
    };
    let new_start = clip(span.start);
    let new_end = clip(span.end);
    if new_end <= new_start {
        None // the span lay entirely within the deletion
    } else {
        Some(Span {
            start: new_start,
            end: new_end,
            strand: span.strand,
        })
    }
}

/// Remaps a [`Location`] for an insertion. Every span is shifted.
fn remap_location_insert(loc: &Location, pos: usize, len: usize) -> Location {
    match loc {
        Location::Single(s) => Location::Single(shift_span_for_insert(*s, pos, len)),
        Location::Join(v) => Location::Join(
            v.iter()
                .map(|s| shift_span_for_insert(*s, pos, len))
                .collect(),
        ),
        Location::Order(v) => Location::Order(
            v.iter()
                .map(|s| shift_span_for_insert(*s, pos, len))
                .collect(),
        ),
        Location::Between { position, strand } => Location::Between {
            // A bond at `position` shifts by `len` if the insertion
            // lies at or before the bond.
            position: if *position >= pos {
                *position + len
            } else {
                *position
            },
            strand: *strand,
        },
    }
}

/// Remaps a [`Location`] for a deletion. Spans fully inside the
/// deletion are dropped; the location is `None` if all spans vanish.
fn remap_location_delete(loc: &Location, del_start: usize, del_end: usize) -> Option<Location> {
    let del_len = del_end - del_start;
    match loc {
        Location::Single(s) => {
            shift_span_for_delete(*s, del_start, del_end).map(Location::Single)
        }
        Location::Join(v) => {
            let kept: Vec<Span> = v
                .iter()
                .filter_map(|s| shift_span_for_delete(*s, del_start, del_end))
                .collect();
            if kept.is_empty() {
                None
            } else if kept.len() == 1 {
                Some(Location::Single(kept[0]))
            } else {
                Some(Location::Join(kept))
            }
        }
        Location::Order(v) => {
            let kept: Vec<Span> = v
                .iter()
                .filter_map(|s| shift_span_for_delete(*s, del_start, del_end))
                .collect();
            if kept.is_empty() {
                None
            } else {
                // Order is semantic; keep it even if only one segment
                // survives, since the surviving entry is still "an
                // ordered list of one".
                Some(Location::Order(kept))
            }
        }
        Location::Between { position, strand } => {
            // A bond inside the deletion is dropped; one outside
            // shifts accordingly.
            if *position > del_start && *position < del_end {
                None
            } else {
                let new_pos = if *position >= del_end {
                    *position - del_len
                } else {
                    *position
                };
                Some(Location::Between {
                    position: new_pos,
                    strand: *strand,
                })
            }
        }
    }
}

/// Inserts `insert` into a [`SeqRecord`] at `pos`, shifting every
/// affected feature's coordinates.
pub fn insert_into_record(rec: &SeqRecord, pos: usize, ins: &Seq) -> Result<SeqRecord> {
    let new_seq = insert(&rec.seq, pos, ins)?;
    let len = ins.len();
    let features: Vec<SeqFeature> = rec
        .features
        .iter()
        .map(|f| SeqFeature {
            feature_type: f.feature_type.clone(),
            location: remap_location_insert(&f.location, pos, len),
            qualifiers: f.qualifiers.clone(),
        })
        .collect();
    Ok(SeqRecord {
        seq: new_seq,
        features,
        ..rec.clone()
    })
}

/// Deletes `[start, end)` from a [`SeqRecord`], clipping or dropping
/// any feature that overlaps the deletion.
pub fn delete_from_record(rec: &SeqRecord, start: usize, end: usize) -> Result<SeqRecord> {
    let new_seq = delete(&rec.seq, start, end)?;
    let features: Vec<SeqFeature> = rec
        .features
        .iter()
        .filter_map(|f| {
            remap_location_delete(&f.location, start, end).map(|loc| SeqFeature {
                feature_type: f.feature_type.clone(),
                location: loc,
                qualifiers: f.qualifiers.clone(),
            })
        })
        .collect();
    Ok(SeqRecord {
        seq: new_seq,
        features,
        ..rec.clone()
    })
}

/// Extracts the sub-record `[start, end)` of a [`SeqRecord`], keeping
/// only the features that fall entirely inside the slice and rebasing
/// their coordinates to the slice origin.
pub fn slice_record(rec: &SeqRecord, start: usize, end: usize) -> Result<SeqRecord> {
    let new_seq = rec.seq.slice(start, end)?;
    let features: Vec<SeqFeature> = rec
        .features
        .iter()
        .filter_map(|f| {
            let loc = &f.location;
            // Keep only features fully within [start, end).
            if loc.min_start() >= start && loc.max_end() <= end {
                let rebase = |s: &Span| Span {
                    start: s.start - start,
                    end: s.end - start,
                    strand: s.strand,
                };
                let new_loc = match loc {
                    Location::Single(s) => Location::Single(rebase(s)),
                    Location::Join(v) => Location::Join(v.iter().map(rebase).collect()),
                    Location::Order(v) => Location::Order(v.iter().map(rebase).collect()),
                    Location::Between { position, strand } => Location::Between {
                        position: position - start,
                        strand: *strand,
                    },
                };
                Some(SeqFeature {
                    feature_type: f.feature_type.clone(),
                    location: new_loc,
                    qualifiers: f.qualifiers.clone(),
                })
            } else {
                None
            }
        })
        .collect();
    Ok(SeqRecord {
        id: rec.id.clone(),
        name: rec.name.clone(),
        description: rec.description.clone(),
        seq: new_seq,
        features,
        annotations: rec.annotations.clone(),
        references: rec.references.clone(),
    })
}

/// Rotates a circular [`SeqRecord`] so that `origin` becomes position
/// 0, wrapping every feature coordinate accordingly.
///
/// Returns [`BioseqError::Invalid`] for a linear record (rotation
/// would split linear features ambiguously).
pub fn rotate_record(rec: &SeqRecord, origin: usize) -> Result<SeqRecord> {
    if !rec.seq.is_circular() {
        return Err(BioseqError::invalid(
            "topology",
            "rotate_record requires a circular record",
        ));
    }
    let n = rec.seq.len();
    if n == 0 {
        return Ok(rec.clone());
    }
    let o = origin % n;
    let new_seq = rec.seq.rotate(o)?;
    let wrap = |c: usize| (c + n - o) % n;
    let features: Vec<SeqFeature> = rec
        .features
        .iter()
        .map(|f| {
            let remap = |s: &Span| {
                let mut ns = wrap(s.start);
                let mut ne = wrap(s.end.saturating_sub(1)) + 1;
                // Guard against a degenerate wrap producing ne <= ns.
                if ne <= ns {
                    std::mem::swap(&mut ns, &mut ne);
                }
                Span {
                    start: ns,
                    end: ne,
                    strand: s.strand,
                }
            };
            let loc = match &f.location {
                Location::Single(s) => Location::Single(remap(s)),
                Location::Join(v) => Location::Join(v.iter().map(remap).collect()),
                Location::Order(v) => Location::Order(v.iter().map(remap).collect()),
                Location::Between { position, strand } => Location::Between {
                    position: (position + n - o) % n,
                    strand: *strand,
                },
            };
            SeqFeature {
                feature_type: f.feature_type.clone(),
                location: loc,
                qualifiers: f.qualifiers.clone(),
            }
        })
        .collect();
    Ok(SeqRecord {
        seq: new_seq,
        features,
        ..rec.clone()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::{SeqKind, Topology};

    fn dna(s: &str) -> Seq {
        Seq::new(SeqKind::Dna, s).unwrap()
    }

    #[test]
    fn insert_basic() {
        let s = dna("AAAATTTT");
        let out = insert(&s, 4, &dna("GG")).unwrap();
        assert_eq!(out.as_str(), "AAAAGGTTTT");
        // Append at the end.
        assert_eq!(insert(&s, 8, &dna("C")).unwrap().as_str(), "AAAATTTTC");
        assert!(insert(&s, 99, &dna("C")).is_err());
    }

    #[test]
    fn delete_basic() {
        let s = dna("AAAATTTT");
        assert_eq!(delete(&s, 2, 6).unwrap().as_str(), "AATT");
        assert!(delete(&s, 6, 2).is_err());
        assert!(delete(&s, 0, 99).is_err());
    }

    #[test]
    fn replace_basic() {
        let s = dna("AAAATTTT");
        // Replace a 2-base stretch with a 4-base one.
        assert_eq!(replace(&s, 4, 6, &dna("GGGG")).unwrap().as_str(), "AAAAGGGGTT");
    }

    #[test]
    fn insert_shifts_downstream_features() {
        let mut rec = SeqRecord::new("r", dna("AAAATTTT"));
        rec.add_feature(SeqFeature::new("upstream", Location::single(0, 2)));
        rec.add_feature(SeqFeature::new("downstream", Location::single(4, 8)));
        let edited = insert_into_record(&rec, 2, &dna("GGG")).unwrap();
        assert_eq!(edited.seq.as_str(), "AAGGGAATTTT");
        // upstream (0..2) is untouched.
        assert_eq!(edited.features[0].location, Location::single(0, 2));
        // downstream (4..8) shifts by 3 -> (7..11).
        assert_eq!(edited.features[1].location, Location::single(7, 11));
    }

    #[test]
    fn delete_clips_overlapping_features() {
        let mut rec = SeqRecord::new("r", dna("AAAATTTTGGGG"));
        rec.add_feature(SeqFeature::new("spanning", Location::single(2, 10)));
        rec.add_feature(SeqFeature::new("inside", Location::single(5, 7)));
        // Delete [4, 8).
        let edited = delete_from_record(&rec, 4, 8).unwrap();
        assert_eq!(edited.seq.as_str(), "AAAAGGGG");
        // The "inside" feature lay entirely within the deletion -> dropped.
        assert_eq!(edited.features.len(), 1);
        // "spanning" (2..10) clips: start 2 stays, end 10 -> 10-4=6.
        assert_eq!(edited.features[0].location, Location::single(2, 6));
    }

    #[test]
    fn slice_record_keeps_contained_features() {
        let mut rec = SeqRecord::new("r", dna("AAAATTTTGGGG"));
        rec.add_feature(SeqFeature::new("contained", Location::single(5, 7)));
        rec.add_feature(SeqFeature::new("straddling", Location::single(2, 9)));
        let sub = slice_record(&rec, 4, 8).unwrap();
        assert_eq!(sub.seq.as_str(), "TTTT");
        // Only "contained" survives; rebased 5..7 -> 1..3.
        assert_eq!(sub.features.len(), 1);
        assert_eq!(sub.features[0].location, Location::single(1, 3));
    }

    #[test]
    fn rotate_record_wraps_features() {
        let mut rec = SeqRecord::new(
            "p",
            Seq::with_topology(SeqKind::Dna, "ATGCAA", Topology::Circular).unwrap(),
        );
        rec.add_feature(SeqFeature::new("f", Location::single(0, 2)));
        let rotated = rotate_record(&rec, 2).unwrap();
        assert_eq!(rotated.seq.as_str(), "GCAAAT");
        // Feature at 0..2 with origin shifted by 2 -> wraps to 4..6.
        assert_eq!(rotated.features[0].location, Location::single(4, 6));
    }

    #[test]
    fn rotate_rejects_linear() {
        let rec = SeqRecord::new("r", dna("ATGC"));
        assert!(rotate_record(&rec, 1).is_err());
    }

    #[test]
    fn kind_mismatch_rejected() {
        let d = dna("ACGT");
        let r = Seq::new(SeqKind::Rna, "ACGU").unwrap();
        assert!(insert(&d, 0, &r).is_err());
        assert!(replace(&d, 0, 1, &r).is_err());
    }
}
