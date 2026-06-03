//! mRNA therapeutic design — Group D (features 16–23).
//!
//! Designing a therapeutic mRNA — the kind made with `m1Ψ` and
//! delivered in a lipid nanoparticle. This block assembles and
//! optimises the five-part construct (cap, 5′UTR, CDS, 3′UTR,
//! poly-A):
//!
//! - [`construct`] — feature 16: the [`construct::MrnaConstruct`] model
//!   and a validating [`construct::MrnaConstructBuilder`].
//! - [`codon`] — feature 17: CDS codon optimisation for expression,
//!   reusing [`valenx_bioseq`]'s codon tables + CAI.
//! - [`utr`] — features 18–19: Kozak-aware 5′UTR design and
//!   AU-rich-element-aware 3′UTR design.
//! - [`structure`] — feature 20: a start-codon-region secondary-
//!   structure check and minimisation, reusing [`valenx_rnastruct`]'s
//!   Zuker MFE folder.
//! - [`uridine`] — feature 21: uridine-content measurement and a
//!   uridine-depleting synonymous-codon pass for pseudouridine
//!   constructs.
//! - [`tailcap`] — feature 22: poly-A tail and 5′ cap-analog
//!   selection.
//! - [`advanced`] — feature 23: self-amplifying-mRNA (saRNA) and
//!   circular-RNA (circRNA) construct-layout design.
//! - [`design`] — the top-level [`design::design_mrna`] driver and the
//!   bundled [`design::MrnaDesignReport`] (the mRNA half of feature 30).
//!
//! ## v1 scope
//!
//! Every score in this block is a transparent heuristic (Kozak match,
//! ARE count, structural openness, end-modification rules), documented
//! as such. The driver chains the passes in a fixed order rather than
//! jointly optimising; each pass's effect is reported. The saRNA /
//! circRNA designers are construct-layout v1s — they order and
//! validate the parts, they do not simulate replication or splicing.

pub mod advanced;
pub mod codon;
pub mod construct;
pub mod design;
pub mod structure;
pub mod tailcap;
pub mod uridine;
pub mod utr;
