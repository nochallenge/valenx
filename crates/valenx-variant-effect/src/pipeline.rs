//! The orchestration pipeline.
//!
//! A [`VariantEffectPipeline`] holds an ordered set of
//! [`VariantEffectPredictor`]s and, given a variant and a reference
//! sequence, builds the [`VariantContext`] once and runs every predictor
//! against it, aggregating the results into a [`VariantReport`].
//!
//! # Predictor-error policy
//!
//! Per-predictor errors are **surfaced, not fatal**. If
//! [`build`](VariantContext::build) succeeds (the variant is valid for the
//! reference), the run always completes: each predictor that returns
//! `Ok` contributes an `EffectPrediction`, and each that returns `Err`
//! contributes a [`PredictionFailure`] (its `Display` text) instead. One
//! broken predictor therefore never hides the verdicts of the others, and
//! [`run`](VariantEffectPipeline::run) returns `Err` only when the
//! variant/reference pair itself is invalid (so there is nothing to
//! predict on).

use valenx_bioseq::seq::Seq;

use crate::error::VariantError;
use crate::predictor::{VariantContext, VariantEffectPredictor};
use crate::report::{PredictionFailure, VariantReport};
use crate::variant::Variant;

/// An ordered pipeline of variant-effect predictors.
#[derive(Default)]
pub struct VariantEffectPipeline {
    predictors: Vec<Box<dyn VariantEffectPredictor>>,
}

impl VariantEffectPipeline {
    /// An empty pipeline (no predictors).
    pub fn new() -> Self {
        VariantEffectPipeline {
            predictors: Vec::new(),
        }
    }

    /// Builder-style: adds a predictor and returns `self`, so predictors
    /// can be chained at construction.
    #[must_use]
    pub fn with_predictor(mut self, p: Box<dyn VariantEffectPredictor>) -> Self {
        self.predictors.push(p);
        self
    }

    /// Adds a predictor to an existing pipeline.
    pub fn add_predictor(&mut self, p: Box<dyn VariantEffectPredictor>) {
        self.predictors.push(p);
    }

    /// The number of predictors in the pipeline.
    pub fn len(&self) -> usize {
        self.predictors.len()
    }

    /// `true` if the pipeline holds no predictors.
    pub fn is_empty(&self) -> bool {
        self.predictors.is_empty()
    }

    /// Runs the pipeline for one variant against one reference sequence.
    ///
    /// Builds the [`VariantContext`] (applying the variant to
    /// `wild_type`), then runs every predictor in order. See the
    /// [module-level error policy](self): per-predictor errors become
    /// [`PredictionFailure`]s in the report; only an invalid
    /// variant/reference pair makes this return `Err`.
    pub fn run(&self, variant: Variant, wild_type: Seq) -> Result<VariantReport, VariantError> {
        let ctx = VariantContext::build(variant, wild_type)?;

        let mut predictions = Vec::new();
        let mut failures = Vec::new();

        for predictor in &self.predictors {
            match predictor.predict(&ctx) {
                Ok(pred) => predictions.push(pred),
                Err(e) => failures.push(PredictionFailure {
                    source: predictor.name().to_string(),
                    error: e.to_string(),
                }),
            }
        }

        Ok(VariantReport::new(ctx.variant, predictions, failures))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::predictor::{EffectCategory, MockPredictor};
    use crate::variant::parse;
    use valenx_bioseq::seq::SeqKind;

    fn dna(s: &str) -> Seq {
        Seq::new(SeqKind::Dna, s).unwrap()
    }

    fn variant() -> Variant {
        parse("c.5G>A").unwrap()
    }

    #[test]
    fn two_predictors_yield_two_predictions() {
        let pipeline = VariantEffectPipeline::new()
            .with_predictor(Box::new(MockPredictor::with_score(
                "Missense",
                EffectCategory::Damaging,
                0.9,
            )))
            .with_predictor(Box::new(MockPredictor::new(
                "Conservation",
                EffectCategory::Damaging,
            )));
        assert_eq!(pipeline.len(), 2);

        let report = pipeline.run(variant(), dna("ATGCGT")).unwrap();
        assert_eq!(report.predictions.len(), 2);
        assert!(report.failures.is_empty());
        assert_eq!(report.consensus, EffectCategory::Damaging);
        // Order is preserved.
        assert_eq!(report.predictions[0].source, "Missense");
        assert_eq!(report.predictions[1].source, "Conservation");
    }

    #[test]
    fn add_predictor_mutates_in_place() {
        let mut pipeline = VariantEffectPipeline::new();
        assert!(pipeline.is_empty());
        pipeline.add_predictor(Box::new(MockPredictor::new("A", EffectCategory::Benign)));
        pipeline.add_predictor(Box::new(MockPredictor::new("B", EffectCategory::Benign)));
        assert_eq!(pipeline.len(), 2);
        let report = pipeline.run(variant(), dna("ATGCGT")).unwrap();
        assert_eq!(report.predictions.len(), 2);
        assert_eq!(report.consensus, EffectCategory::Benign);
    }

    #[test]
    fn erroring_predictor_is_surfaced_not_fatal() {
        // One good predictor, one that always errors. The run must succeed,
        // returning the good prediction AND recording the failure.
        let pipeline = VariantEffectPipeline::new()
            .with_predictor(Box::new(MockPredictor::new(
                "Good",
                EffectCategory::Damaging,
            )))
            .with_predictor(Box::new(MockPredictor::failing(
                "Broken",
                VariantError::Bioseq("model unavailable".into()),
            )));

        let report = pipeline.run(variant(), dna("ATGCGT")).unwrap();
        assert_eq!(report.predictions.len(), 1);
        assert_eq!(report.predictions[0].source, "Good");
        assert_eq!(report.failures.len(), 1);
        assert_eq!(report.failures[0].source, "Broken");
        assert!(
            report.failures[0].error.contains("model unavailable"),
            "got: {}",
            report.failures[0].error
        );
        // Consensus is computed over the surviving prediction only.
        assert_eq!(report.consensus, EffectCategory::Damaging);
    }

    #[test]
    fn all_predictors_erroring_yields_empty_predictions_uncertain_consensus() {
        let pipeline = VariantEffectPipeline::new()
            .with_predictor(Box::new(MockPredictor::failing(
                "B1",
                VariantError::Bioseq("x".into()),
            )))
            .with_predictor(Box::new(MockPredictor::failing(
                "B2",
                VariantError::Bioseq("y".into()),
            )));
        let report = pipeline.run(variant(), dna("ATGCGT")).unwrap();
        assert!(report.predictions.is_empty());
        assert_eq!(report.failures.len(), 2);
        assert_eq!(report.consensus, EffectCategory::Uncertain);
    }

    #[test]
    fn empty_pipeline_yields_report_with_no_predictions() {
        let pipeline = VariantEffectPipeline::new();
        let report = pipeline.run(variant(), dna("ATGCGT")).unwrap();
        assert!(report.predictions.is_empty());
        assert!(report.failures.is_empty());
        assert_eq!(report.consensus, EffectCategory::Uncertain);
        assert_eq!(report.variant, variant());
    }

    #[test]
    fn invalid_variant_reference_pair_aborts_the_run() {
        // Wild-type mismatch (position 4 is C, not G) means there is
        // nothing to predict on -> the run returns Err, not a report.
        let pipeline = VariantEffectPipeline::new().with_predictor(Box::new(MockPredictor::new(
            "Good",
            EffectCategory::Damaging,
        )));
        let result = pipeline.run(parse("c.4G>A").unwrap(), dna("ATGCGT"));
        assert!(matches!(result, Err(VariantError::WildTypeMismatch { .. })));
    }

    #[test]
    fn default_is_empty_pipeline() {
        let pipeline = VariantEffectPipeline::default();
        assert!(pipeline.is_empty());
    }
}
