//! The pluggable predictor seam.
//!
//! A [`VariantEffectPredictor`] takes a fully-built [`VariantContext`]
//! (the variant plus its wild-type and mutant sequences) and returns an
//! [`EffectPrediction`]. Real predictors wrap the existing tool adapters
//! (AlphaMissense and friends); the bundled [`MockPredictor`] returns a
//! configured prediction deterministically so the pipeline and report
//! logic can be tested without any external tool.

use serde::Serialize;

use valenx_bioseq::seq::Seq;

use crate::apply::apply;
use crate::error::VariantError;
use crate::variant::Variant;

/// Everything a predictor needs about one variant: the parsed variant and
/// both the wild-type and mutant sequences.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VariantContext {
    /// The variant being assessed.
    pub variant: Variant,
    /// The reference (wild-type) sequence the variant was applied to.
    pub wild_type: Seq,
    /// The mutant sequence produced by applying the variant.
    pub mutant: Seq,
}

impl VariantContext {
    /// Builds a context by applying `variant` to `wild_type`.
    ///
    /// Returns whatever [`apply`] returns on a kind / range / wild-type
    /// mismatch, so an invalid variant/reference pair is rejected up front
    /// rather than reaching the predictors.
    pub fn build(variant: Variant, wild_type: Seq) -> Result<Self, VariantError> {
        let mutant = apply(&variant, &wild_type)?;
        Ok(VariantContext {
            variant,
            wild_type,
            mutant,
        })
    }
}

/// The qualitative effect call.
///
/// Ordered from least to most damaging so the discriminant can double as
/// a severity rank; the explicit values keep the wire format stable.
///
/// Marked `#[non_exhaustive]`: this crate is an extensible orchestration
/// layer and more effect categories may be added in future, so downstream
/// matches must include a wildcard arm.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize)]
#[non_exhaustive]
pub enum EffectCategory {
    /// Predicted not to affect function.
    Benign = 0,
    /// Insufficient or conflicting evidence to call.
    Uncertain = 1,
    /// Predicted to damage / disrupt function (the "likely pathogenic"
    /// end of the scale).
    Damaging = 2,
}

impl EffectCategory {
    /// A short lowercase label (`"benign"`, `"uncertain"`, `"damaging"`).
    pub fn label(self) -> &'static str {
        match self {
            EffectCategory::Benign => "benign",
            EffectCategory::Uncertain => "uncertain",
            EffectCategory::Damaging => "damaging",
        }
    }
}

/// One predictor's verdict on a variant.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct EffectPrediction {
    /// Name of the predictor that produced this (its
    /// [`VariantEffectPredictor::name`]).
    pub source: String,
    /// Optional numeric score, in whatever scale the predictor uses
    /// (`None` if the predictor is purely categorical).
    pub score: Option<f64>,
    /// The qualitative call.
    pub category: EffectCategory,
    /// Free-text notes (model version, caveats, raw label, …).
    pub notes: String,
}

/// A pluggable variant-effect predictor.
///
/// Implementors must be `Send + Sync` so a [`VariantEffectPipeline`] can
/// hold them as `Box<dyn VariantEffectPredictor>` and (in future) run them
/// concurrently.
///
/// [`VariantEffectPipeline`]: crate::pipeline::VariantEffectPipeline
pub trait VariantEffectPredictor: Send + Sync {
    /// A stable, human-readable name for this predictor (used as the
    /// [`EffectPrediction::source`]).
    fn name(&self) -> &str;

    /// Predicts the effect of the variant described by `ctx`.
    fn predict(&self, ctx: &VariantContext) -> Result<EffectPrediction, VariantError>;
}

/// How a [`MockPredictor`] behaves when [`predict`](VariantEffectPredictor::predict)
/// is called.
#[derive(Clone, Debug, PartialEq)]
pub enum MockBehavior {
    /// Always return a prediction with this category and score.
    Fixed {
        /// Category the mock reports.
        category: EffectCategory,
        /// Score the mock reports.
        score: Option<f64>,
    },
    /// Always fail with a [`VariantError`] (lets pipeline tests exercise
    /// the predictor-error path deterministically).
    Error(VariantError),
}

/// A deterministic predictor for tests and examples.
///
/// Constructed with a [`MockBehavior`], it returns the configured
/// prediction (or error) for every context, ignoring the context's
/// contents. This is the reference implementation of
/// [`VariantEffectPredictor`].
#[derive(Clone, Debug)]
pub struct MockPredictor {
    name: String,
    behavior: MockBehavior,
}

impl MockPredictor {
    /// A mock that always returns `category` with no numeric score.
    pub fn new(name: impl Into<String>, category: EffectCategory) -> Self {
        MockPredictor {
            name: name.into(),
            behavior: MockBehavior::Fixed {
                category,
                score: None,
            },
        }
    }

    /// A mock that always returns `category` with the given `score`.
    pub fn with_score(
        name: impl Into<String>,
        category: EffectCategory,
        score: f64,
    ) -> Self {
        MockPredictor {
            name: name.into(),
            behavior: MockBehavior::Fixed {
                category,
                score: Some(score),
            },
        }
    }

    /// A mock that always fails with `error` — for exercising the
    /// pipeline's predictor-error policy.
    pub fn failing(name: impl Into<String>, error: VariantError) -> Self {
        MockPredictor {
            name: name.into(),
            behavior: MockBehavior::Error(error),
        }
    }
}

impl VariantEffectPredictor for MockPredictor {
    fn name(&self) -> &str {
        &self.name
    }

    fn predict(&self, _ctx: &VariantContext) -> Result<EffectPrediction, VariantError> {
        match &self.behavior {
            MockBehavior::Fixed { category, score } => Ok(EffectPrediction {
                source: self.name.clone(),
                score: *score,
                category: *category,
                notes: format!("mock predictor `{}`", self.name),
            }),
            MockBehavior::Error(e) => Err(e.clone()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::parse;
    use valenx_bioseq::seq::SeqKind;

    fn dna(s: &str) -> Seq {
        Seq::new(SeqKind::Dna, s).unwrap()
    }

    fn ctx() -> VariantContext {
        VariantContext::build(parse("c.5G>A").unwrap(), dna("ATGCGT")).unwrap()
    }

    #[test]
    fn context_build_populates_wild_type_and_mutant() {
        let c = ctx();
        assert_eq!(c.wild_type.as_str(), "ATGCGT");
        assert_eq!(c.mutant.as_str(), "ATGCAT");
        assert_eq!(c.variant, parse("c.5G>A").unwrap());
    }

    #[test]
    fn context_build_propagates_apply_errors() {
        // Wild-type mismatch: position 4 is C, not G.
        let err = VariantContext::build(parse("c.4G>A").unwrap(), dna("ATGCGT"));
        assert!(matches!(err, Err(VariantError::WildTypeMismatch { .. })));
    }

    #[test]
    fn effect_category_ordering_and_labels() {
        assert!(EffectCategory::Benign < EffectCategory::Uncertain);
        assert!(EffectCategory::Uncertain < EffectCategory::Damaging);
        assert_eq!(EffectCategory::Damaging.label(), "damaging");
    }

    #[test]
    fn mock_returns_configured_prediction() {
        let p = MockPredictor::with_score("MockMissense", EffectCategory::Damaging, 0.97);
        let pred = p.predict(&ctx()).unwrap();
        assert_eq!(pred.source, "MockMissense");
        assert_eq!(pred.category, EffectCategory::Damaging);
        assert_eq!(pred.score, Some(0.97));
        assert!(pred.notes.contains("MockMissense"));
    }

    #[test]
    fn mock_without_score_has_none() {
        let p = MockPredictor::new("Cat", EffectCategory::Benign);
        let pred = p.predict(&ctx()).unwrap();
        assert_eq!(pred.score, None);
        assert_eq!(pred.category, EffectCategory::Benign);
    }

    #[test]
    fn mock_failing_returns_error() {
        let p = MockPredictor::failing(
            "Broken",
            VariantError::Bioseq("model unavailable".into()),
        );
        assert!(matches!(p.predict(&ctx()), Err(VariantError::Bioseq(_))));
    }

    #[test]
    fn usable_as_trait_object() {
        let preds: Vec<Box<dyn VariantEffectPredictor>> = vec![
            Box::new(MockPredictor::new("A", EffectCategory::Benign)),
            Box::new(MockPredictor::new("B", EffectCategory::Damaging)),
        ];
        let c = ctx();
        let cats: Vec<_> = preds
            .iter()
            .map(|p| p.predict(&c).unwrap().category)
            .collect();
        assert_eq!(cats, vec![EffectCategory::Benign, EffectCategory::Damaging]);
        // names accessible through the trait object
        assert_eq!(preds[0].name(), "A");
    }
}
