//! # valenx-binder-score
//!
//! Fold a candidate binder's heterogeneous evidence — a **binding-affinity
//! estimate**, a **developability** summary, and the **worst safety severity** —
//! into one comparable **binder-quality score** in `[0, 1]`, with every
//! component kept visible. The single number a pipeline ranks on after scoring,
//! developability, and safety screening have each run.
//!
//! ## What
//!
//! - [`BinderInputs`] — the (all-optional) evidence channels.
//! - [`BinderWeights`] — per-channel weights (default `1.0`; `0` drops a channel).
//! - [`score()`] / [`BinderScore`] — the fused `[0, 1]` value plus the
//!   normalized, weighted [`BinderScore::components`].
//!
//! ## Model
//!
//! Each present channel maps to `[0, 1]` (higher = better): binding ΔG through a
//! monotone "more-negative-is-better" transform; developability passed through
//! (already `[0, 1]`); safety severity `0..=4` as `1 − severity/4` (0 = no
//! flags → 1, 4 = critical → 0). The score is the weighted mean of the present
//! channels. It is monotone — better binding or lower severity never lowers the
//! score.
//!
//! ## The rule it keeps
//!
//! A high score is **not** a safety clearance. [`BinderScore::requires_review`]
//! is always `true`: this ranks candidates for a human to triage, it never
//! approves one for the wet lab.
//!
//! ## Honest scope
//!
//! Research/educational. This is a transparent **selection heuristic** that
//! combines upstream estimates — it is only as good as its inputs and the
//! (arbitrary, documented) normalization, and it makes no validated claim about
//! a candidate's real quality, affinity, or safety. Calibrate the ΔG transform
//! against data before reading the score as anything but a rank.
//!
//! ## Example
//!
//! ```
//! use valenx_binder_score::{score, BinderInputs, BinderWeights};
//!
//! let strong = BinderInputs { dg_bind_kcal: Some(-12.0), developability: Some(0.8), safety_severity: Some(0) };
//! let s = score(&strong, &BinderWeights::default()).unwrap();
//! assert!(s.value > 0.6 && s.requires_review);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod score;

pub use error::BinderError;
pub use score::{score, BinderInputs, BinderScore, BinderWeights};
