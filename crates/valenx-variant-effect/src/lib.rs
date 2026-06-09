//! # valenx-variant-effect — variant-effect orchestration
//!
//! The integration layer that turns *a genetic variant* into *a
//! predicted-effect report*. It is intentionally thin: parsing,
//! sequence-application and aggregation live here, while the actual
//! effect predictions are delegated to pluggable
//! [`VariantEffectPredictor`]
//! implementations (the real ones are backed by the existing tool
//! adapters; tests use [`MockPredictor`]).
//!
//! This crate is a **pure library**: it writes no files. The
//! [`run`](pipeline::VariantEffectPipeline::run) entry point returns a
//! [`VariantReport`]; persisting it (e.g. as JSON)
//! is the caller's responsibility.
//!
//! ## Pipeline
//!
//! The flow is: parse a [`Variant`] → build a
//! [`VariantEffectPipeline`] of [`VariantEffectPredictor`]s →
//! [`run`](VariantEffectPipeline::run) it against a reference
//! [`Seq`](valenx_bioseq::Seq) → read the [`VariantReport`]. The pipeline
//! applies the variant to the reference once (deriving the mutant
//! sequence), runs each predictor, and aggregates a consensus. A
//! predictor that errors is surfaced as a
//! [`PredictionFailure`] in the report rather than aborting the run.
//!
//! ## Example
//!
//! ```rust
//! use valenx_bioseq::{Seq, SeqKind};
//! use valenx_variant_effect::{
//!     variant, EffectCategory, MockPredictor, VariantEffectPipeline,
//! };
//!
//! // 1. Parse a coding substitution (c.5G>A: CGT -> CAT, Arg -> His).
//! let v = variant::parse("c.5G>A").unwrap();
//!
//! // 2. Build a pipeline. Real predictors wrap the tool adapters; here
//! //    we plug in deterministic mocks via the builder.
//! let pipeline = VariantEffectPipeline::new()
//!     .with_predictor(Box::new(MockPredictor::with_score(
//!         "Missense",
//!         EffectCategory::Damaging,
//!         0.97,
//!     )))
//!     .with_predictor(Box::new(MockPredictor::new(
//!         "Conservation",
//!         EffectCategory::Damaging,
//!     )));
//!
//! // 3. Run it against a reference coding sequence (ATG CGT = Met Arg).
//! let reference = Seq::new(SeqKind::Dna, "ATGCGT").unwrap();
//! let report = pipeline.run(v, reference).unwrap();
//!
//! // 4. Inspect the aggregated report.
//! assert_eq!(report.predictions.len(), 2);
//! assert_eq!(report.consensus, EffectCategory::Damaging);
//! println!("{}", report.summary()); // human-readable; no file is written
//! ```

#![forbid(unsafe_code)]
// `missing_docs` is enforced via `[lints.rust]` in Cargo.toml (matching
// the sibling crates' convention); no redundant crate-level attr here.

pub mod apply;
pub mod error;
pub mod pipeline;
pub mod predictor;
pub mod report;
pub mod variant;

pub use apply::{apply, translate_coding};
pub use error::VariantError;
pub use pipeline::VariantEffectPipeline;
pub use predictor::{
    EffectCategory, EffectPrediction, MockPredictor, VariantContext, VariantEffectPredictor,
};
pub use report::{consensus, PredictionFailure, VariantReport};
pub use variant::{parse, transition_transversion_ratio, AminoAcid, Nucleotide, Variant};

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, VariantError>;
