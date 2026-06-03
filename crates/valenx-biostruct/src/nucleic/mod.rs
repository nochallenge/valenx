//! Nucleic-acid geometry: base-pair detection, 3DNA-class base-pair
//! and step parameters, groove widths and helical-axis fitting.
//!
//! - [`basepair`] — Watson-Crick / wobble base-pair detection by
//!   hydrogen-bond geometry.
//! - [`params`] — base reference frames and the six base-pair
//!   parameters / six step parameters.
//! - [`grooves`] — major- and minor-groove widths.
//! - [`helix`] — global straight axis ([`fit_helical_axis`]) and a
//!   Curves+-class **curved** spline axis with curvature analysis
//!   ([`fit_curved_axis`]).

pub mod basepair;
pub mod grooves;
pub mod helix;
pub mod params;

pub use basepair::{detect_base_pairs, BasePair, BasePairKind};
pub use grooves::{groove_widths, GrooveProfile, GrooveWidth};
pub use helix::{
    fit_curved_axis, fit_helical_axis, ideal_bent_helix_frames, ideal_helix_centers,
    ideal_helix_frames, CurvedHelicalAxis, HelicalAxis, SplineSegment,
};
pub use params::{
    base_frame, base_pair_parameters, pair_mid_frame, step_parameters, BaseFrame,
    BasePairParameters, StepParameters,
};
