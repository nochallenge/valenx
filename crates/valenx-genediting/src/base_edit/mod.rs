//! Base editing — Group B (features 6–10).
//!
//! Base editors install a single point mutation without a
//! double-strand break — a deaminase fused to a Cas nickase. This
//! block designs base-editing experiments:
//!
//! - [`editor`] — feature 6: the base-editor database (CBE: BE3,
//!   BE4max, AncBE4max; ABE: ABE7.10, ABE8e) with PAM, activity window
//!   and target transition.
//! - [`design`] — features 7–10: the editable-site finder
//!   ([`design::find_editable_sites`]), SNV-correcting guide design
//!   ([`design::design_base_edit`]), editing-window / bystander
//!   analysis ([`design::analyze_window`]) and the product-purity
//!   heuristic ([`design::product_purity`]).
//!
//! ## v1 scope
//!
//! The product-purity score is a transparent feature-weighted
//! heuristic (bystander count, window centring, CBE `GpC` context),
//! documented as such — not a trained model. Activity windows are the
//! literature consensus; real activity tails are sequence-context
//! dependent.

pub mod design;
pub mod editor;
