//! The aggregated variant-effect report.
//!
//! A [`VariantReport`] bundles the variant, every predictor's
//! [`EffectPrediction`], any predictor [`failures`](VariantReport::failures),
//! and a [`consensus`](VariantReport::consensus) call derived from the
//! successful predictions. [`VariantReport::summary`] renders it as a
//! human-readable string — it returns the text; it never writes a file.

use serde::Serialize;

use crate::predictor::{EffectCategory, EffectPrediction};
use crate::variant::Variant;

/// A predictor that failed to produce a prediction.
///
/// The pipeline surfaces failures here instead of aborting the whole run,
/// so a single broken predictor never hides the others' verdicts.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct PredictionFailure {
    /// Name of the predictor that failed.
    pub source: String,
    /// The error's `Display` text.
    pub error: String,
}

/// Computes the consensus category over a set of predictions.
///
/// Rule: the most frequent [`EffectCategory`] wins. A tie for the top
/// count — including the empty-input case (no predictions at all) —
/// resolves to [`EffectCategory::Uncertain`], the conservative middle
/// call. A single prediction trivially yields its own category.
pub fn consensus(predictions: &[EffectPrediction]) -> EffectCategory {
    let mut counts = [0usize; 3]; // [Benign, Uncertain, Damaging]
    for p in predictions {
        counts[p.category as usize] += 1;
    }
    let max = counts.iter().copied().max().unwrap_or(0);
    if max == 0 {
        return EffectCategory::Uncertain; // no predictions
    }
    let leaders = counts.iter().filter(|&&c| c == max).count();
    if leaders > 1 {
        return EffectCategory::Uncertain; // tie -> conservative
    }
    // Exactly one leader; map its index back to the category.
    match counts.iter().position(|&c| c == max).unwrap() {
        0 => EffectCategory::Benign,
        1 => EffectCategory::Uncertain,
        _ => EffectCategory::Damaging,
    }
}

/// The aggregated report for one variant across all predictors.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct VariantReport {
    /// The variant assessed.
    pub variant: Variant,
    /// Each predictor's successful prediction, in the order the pipeline
    /// ran them.
    pub predictions: Vec<EffectPrediction>,
    /// Predictors that errored (surfaced, not fatal). Empty on a fully
    /// successful run.
    pub failures: Vec<PredictionFailure>,
    /// The consensus call over [`predictions`](Self::predictions).
    pub consensus: EffectCategory,
}

impl VariantReport {
    /// Builds a report from a variant and its predictions/failures,
    /// computing the consensus over the successful predictions.
    pub fn new(
        variant: Variant,
        predictions: Vec<EffectPrediction>,
        failures: Vec<PredictionFailure>,
    ) -> Self {
        let consensus = consensus(&predictions);
        VariantReport {
            variant,
            predictions,
            failures,
            consensus,
        }
    }

    /// A human-readable multi-line summary of the report.
    ///
    /// Returns the rendered string; this function performs no I/O. The
    /// caller decides whether to print it, log it, or persist it.
    pub fn summary(&self) -> String {
        use std::fmt::Write as _;
        let mut s = String::new();
        let _ = writeln!(s, "Variant: {}", describe_variant(&self.variant));
        let _ = writeln!(s, "Consensus: {}", self.consensus.label());
        if self.predictions.is_empty() {
            let _ = writeln!(s, "Predictions: (none)");
        } else {
            let _ = writeln!(s, "Predictions ({}):", self.predictions.len());
            for p in &self.predictions {
                match p.score {
                    Some(score) => {
                        let _ = writeln!(
                            s,
                            "  - {}: {} (score {:.3}) — {}",
                            p.source,
                            p.category.label(),
                            score,
                            p.notes
                        );
                    }
                    None => {
                        let _ =
                            writeln!(s, "  - {}: {} — {}", p.source, p.category.label(), p.notes);
                    }
                }
            }
        }
        if !self.failures.is_empty() {
            let _ = writeln!(s, "Failures ({}):", self.failures.len());
            for f in &self.failures {
                let _ = writeln!(s, "  - {}: {}", f.source, f.error);
            }
        }
        s
    }
}

/// Renders a [`Variant`] back into its HGVS-style string.
fn describe_variant(v: &Variant) -> String {
    match v {
        Variant::ProteinSub { wt, pos, mt } => {
            format!("p.{}{}{}", wt.as_char(), pos, mt.as_char())
        }
        Variant::CodingSub { pos, wt, mt } => {
            format!("c.{}{}>{}", pos, wt.as_char(), mt.as_char())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::parse;

    fn pred(source: &str, category: EffectCategory) -> EffectPrediction {
        EffectPrediction {
            source: source.to_string(),
            score: None,
            category,
            notes: String::new(),
        }
    }

    // ---- consensus rule --------------------------------------------------

    #[test]
    fn consensus_empty_is_uncertain() {
        assert_eq!(consensus(&[]), EffectCategory::Uncertain);
    }

    #[test]
    fn consensus_single_is_that_category() {
        assert_eq!(
            consensus(&[pred("a", EffectCategory::Damaging)]),
            EffectCategory::Damaging
        );
        assert_eq!(
            consensus(&[pred("a", EffectCategory::Benign)]),
            EffectCategory::Benign
        );
    }

    #[test]
    fn consensus_majority_wins() {
        let v = [
            pred("a", EffectCategory::Damaging),
            pred("b", EffectCategory::Damaging),
            pred("c", EffectCategory::Benign),
        ];
        assert_eq!(consensus(&v), EffectCategory::Damaging);
    }

    #[test]
    fn consensus_tie_is_uncertain() {
        // 1 Benign vs 1 Damaging -> tie -> Uncertain.
        let v = [
            pred("a", EffectCategory::Benign),
            pred("b", EffectCategory::Damaging),
        ];
        assert_eq!(consensus(&v), EffectCategory::Uncertain);
    }

    #[test]
    fn consensus_three_way_tie_is_uncertain() {
        let v = [
            pred("a", EffectCategory::Benign),
            pred("b", EffectCategory::Uncertain),
            pred("c", EffectCategory::Damaging),
        ];
        assert_eq!(consensus(&v), EffectCategory::Uncertain);
    }

    #[test]
    fn consensus_uncertain_can_win_outright() {
        let v = [
            pred("a", EffectCategory::Uncertain),
            pred("b", EffectCategory::Uncertain),
            pred("c", EffectCategory::Benign),
        ];
        assert_eq!(consensus(&v), EffectCategory::Uncertain);
    }

    #[test]
    fn consensus_two_two_one_split_ties_to_uncertain() {
        // 2 Benign + 2 Damaging + 1 Uncertain: the top count (2) is shared
        // by two categories, so there is no single leader -> tie ->
        // Uncertain (the conservative call), even though Uncertain itself
        // is the minority here.
        let v = [
            pred("a", EffectCategory::Benign),
            pred("b", EffectCategory::Benign),
            pred("c", EffectCategory::Damaging),
            pred("d", EffectCategory::Damaging),
            pred("e", EffectCategory::Uncertain),
        ];
        assert_eq!(consensus(&v), EffectCategory::Uncertain);
    }

    // ---- report construction + summary ----------------------------------

    #[test]
    fn report_new_computes_consensus() {
        let v = parse("p.R273H").unwrap();
        let report = VariantReport::new(
            v,
            vec![
                pred("a", EffectCategory::Damaging),
                pred("b", EffectCategory::Damaging),
            ],
            vec![],
        );
        assert_eq!(report.consensus, EffectCategory::Damaging);
        assert_eq!(report.predictions.len(), 2);
        assert!(report.failures.is_empty());
    }

    #[test]
    fn summary_mentions_variant_consensus_and_each_predictor() {
        let v = parse("p.R273H").unwrap();
        let report = VariantReport::new(
            v,
            vec![
                EffectPrediction {
                    source: "MockMissense".into(),
                    score: Some(0.97),
                    category: EffectCategory::Damaging,
                    notes: "demo".into(),
                },
                pred("Conservation", EffectCategory::Benign),
            ],
            vec![PredictionFailure {
                source: "Broken".into(),
                error: "model unavailable".into(),
            }],
        );
        let text = report.summary();
        assert!(text.contains("p.R273H"), "got:\n{text}");
        assert!(text.contains("uncertain"), "tie consensus; got:\n{text}");
        assert!(text.contains("MockMissense"), "got:\n{text}");
        assert!(text.contains("0.97"), "score rendered; got:\n{text}");
        assert!(text.contains("Conservation"), "got:\n{text}");
        assert!(text.contains("Broken"), "failure surfaced; got:\n{text}");
        assert!(text.contains("model unavailable"), "got:\n{text}");
    }

    #[test]
    fn summary_handles_empty_predictions() {
        let v = parse("c.817C>T").unwrap();
        let report = VariantReport::new(v, vec![], vec![]);
        let text = report.summary();
        assert!(text.contains("c.817C>T"), "got:\n{text}");
        assert!(text.contains("(none)"), "got:\n{text}");
        assert!(text.contains("uncertain"), "got:\n{text}");
    }

    // ---- describe_variant round-trip ------------------------------------

    #[test]
    fn describe_variant_round_trips_through_parse() {
        // 1-letter protein and coding forms render to a string that parses
        // back to the identical variant.
        for input in ["p.R273H", "c.817C>T", "p.A100A"] {
            let v = parse(input).unwrap();
            let rendered = describe_variant(&v);
            assert_eq!(rendered, input, "render must match 1-letter input");
            assert_eq!(parse(&rendered).unwrap(), v, "rendered form must re-parse");
        }
    }

    #[test]
    fn describe_variant_normalizes_three_letter_to_one_letter() {
        // 3-letter input is normalized to 1-letter internally, so the
        // rendered form is the 1-letter variant (which still re-parses to
        // the same Variant).
        let v = parse("p.Arg273His").unwrap();
        let rendered = describe_variant(&v);
        assert_eq!(rendered, "p.R273H", "3-letter input normalizes to 1-letter");
        assert_eq!(parse(&rendered).unwrap(), v);
    }

    #[test]
    fn report_serializes_to_json() {
        let v = parse("p.R273H").unwrap();
        let report = VariantReport::new(
            v,
            vec![EffectPrediction {
                source: "M".into(),
                score: Some(0.5),
                category: EffectCategory::Damaging,
                notes: "n".into(),
            }],
            vec![],
        );
        // serde derive must produce valid JSON via a serde data format.
        let json = serde_json::to_string(&report).expect("report serializes");
        assert!(json.contains("\"consensus\""), "got: {json}");
        assert!(json.contains("Damaging"), "got: {json}");
        assert!(json.contains("ProteinSub"), "variant tag; got: {json}");
    }
}
