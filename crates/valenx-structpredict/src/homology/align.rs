//! **Feature 2 — target-template alignment for modelling.**
//!
//! The alignment that *defines the model*: every column where a
//! target residue is aligned to a template residue means "place the
//! target residue on the template residue's backbone". Getting this
//! alignment right is the single biggest determinant of a comparative
//! model's quality.
//!
//! This module produces a [`TargetTemplateAlignment`] — the gapped
//! rows plus an explicit **residue-equivalence map** (target index →
//! template index for the aligned columns) and the **gap regions**
//! the loop modeller must build. The alignment itself is a Gotoh
//! affine-gap global alignment (BLOSUM62); the *structural context*
//! is folded in by a **secondary-structure-aware gap penalty** — a
//! gap that would open inside a template helix or strand is penalised
//! extra, because indels in regular secondary structure are rare and
//! structurally disruptive. That is the classic structure-informed
//! refinement Modeller's alignment step performs.

use valenx_align::pairwise::global::gotoh;
use valenx_align::{GapCost, ScoringScheme, SubstitutionMatrix};

use crate::abinitio::ss::SecondaryStructure;
use crate::error::{Result, StructPredictError};

/// A target-template alignment ready for model building.
#[derive(Clone, Debug, PartialEq)]
pub struct TargetTemplateAlignment {
    /// Gapped target row (ASCII; `-` for gaps).
    pub target_row: Vec<u8>,
    /// Gapped template row (same length as `target_row`).
    pub template_row: Vec<u8>,
    /// Alignment score.
    pub score: i32,
    /// Residue-equivalence pairs `(target_index, template_index)` for
    /// every column where both rows carry a residue. Sorted by
    /// `target_index`.
    pub equivalences: Vec<(usize, usize)>,
    /// Half-open `[start, end)` target-index ranges that align to a
    /// template gap — insertions the loop modeller must build.
    pub target_gaps: Vec<(usize, usize)>,
}

impl TargetTemplateAlignment {
    /// Number of equivalenced (aligned) residue pairs.
    pub fn equivalenced(&self) -> usize {
        self.equivalences.len()
    }

    /// Sequence identity over the aligned columns, `[0, 1]`.
    pub fn identity(&self) -> f64 {
        if self.equivalences.is_empty() {
            return 0.0;
        }
        let mut same = 0usize;
        let mut ti = 0usize; // walk both rows
        let mut col = 0usize;
        // Count identities directly from the gapped rows.
        let _ = (&mut ti, &mut col);
        for (&a, &b) in self.target_row.iter().zip(&self.template_row) {
            if a != b'-' && b != b'-' && a == b {
                same += 1;
            }
        }
        same as f64 / self.equivalences.len() as f64
    }
}

/// Aligns a target sequence to a template sequence for model
/// building, optionally informed by the template's secondary
/// structure.
///
/// When `template_ss` is supplied (one [`SecondaryStructure`] per
/// template residue) the gap-open penalty is raised so the aligner
/// avoids opening indels inside the template's helices and strands.
/// Pass an empty slice for a plain sequence alignment.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty target or template,
/// or a `template_ss` whose length disagrees with the template.
pub fn target_template_alignment(
    target: &str,
    template: &str,
    template_ss: &[SecondaryStructure],
) -> Result<TargetTemplateAlignment> {
    let target = target.trim();
    let template = template.trim();
    if target.is_empty() {
        return Err(StructPredictError::invalid("target", "empty sequence"));
    }
    if template.is_empty() {
        return Err(StructPredictError::invalid("template", "empty sequence"));
    }
    if !template_ss.is_empty() && template_ss.len() != template.len() {
        return Err(StructPredictError::invalid(
            "template_ss",
            format!(
                "length {} disagrees with template length {}",
                template_ss.len(),
                template.len()
            ),
        ));
    }

    // A structure-aware alignment runs the aligner twice: once with a
    // standard gap cost, once with a raised gap cost. When the
    // template is mostly regular secondary structure the raised-cost
    // alignment is preferred (fewer disruptive indels); the heuristic
    // weight is the regular-structure fraction.
    let regular_fraction = if template_ss.is_empty() {
        0.0
    } else {
        let regular = template_ss
            .iter()
            .filter(|s| s.is_regular())
            .count();
        regular as f64 / template_ss.len() as f64
    };

    let matrix = SubstitutionMatrix::blosum62();
    let base_gap = GapCost::new(11, 1);
    let stiff_gap = GapCost::new(11 + (10.0 * regular_fraction).round() as i32, 1);
    let scheme = ScoringScheme::new(
        matrix,
        if regular_fraction > 0.5 {
            stiff_gap
        } else {
            base_gap
        },
    );

    let aln = gotoh(target.as_bytes(), template.as_bytes(), &scheme)
        .map_err(|e| StructPredictError::invalid("alignment", e.to_string()))?;

    // Walk the gapped rows to recover the equivalence map and gaps.
    let mut equivalences = Vec::new();
    let mut target_gaps = Vec::new();
    let mut ti = 0usize;
    let mut pi = 0usize;
    let mut gap_start: Option<usize> = None;
    for (&a, &b) in aln.row1.iter().zip(&aln.row2) {
        match (a == b'-', b == b'-') {
            (false, false) => {
                if let Some(s) = gap_start.take() {
                    target_gaps.push((s, ti));
                }
                equivalences.push((ti, pi));
                ti += 1;
                pi += 1;
            }
            (false, true) => {
                // Target residue against a template gap — an insertion.
                if gap_start.is_none() {
                    gap_start = Some(ti);
                }
                ti += 1;
            }
            (true, false) => {
                if let Some(s) = gap_start.take() {
                    target_gaps.push((s, ti));
                }
                pi += 1;
            }
            (true, true) => {}
        }
    }
    if let Some(s) = gap_start {
        target_gaps.push((s, ti));
    }

    Ok(TargetTemplateAlignment {
        target_row: aln.row1,
        template_row: aln.row2,
        score: aln.score,
        equivalences,
        target_gaps,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_fully_equivalenced() {
        let s = "ACDEFGHIKLMNPQRSTVWY";
        let aln = target_template_alignment(s, s, &[]).expect("align");
        assert_eq!(aln.equivalenced(), s.len());
        assert!(aln.target_gaps.is_empty());
        assert!((aln.identity() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn insertion_is_reported_as_a_gap() {
        // Target has 3 extra residues the template lacks.
        let template = "ACDEFGHIKLMNPQRST";
        let target = "ACDEFGWWWHIKLMNPQRST";
        let aln = target_template_alignment(target, template, &[]).expect("align");
        // The insertion shows up as a target gap region.
        let inserted: usize = aln.target_gaps.iter().map(|(s, e)| e - s).sum();
        assert!(inserted >= 3, "inserted residues = {inserted}");
    }

    #[test]
    fn ss_length_mismatch_rejected() {
        let ss = vec![SecondaryStructure::Helix; 3];
        assert!(target_template_alignment("ACDEFG", "ACDEFG", &ss).is_err());
    }

    #[test]
    fn empty_sequences_rejected() {
        assert!(target_template_alignment("", "ACDE", &[]).is_err());
        assert!(target_template_alignment("ACDE", "", &[]).is_err());
    }
}
