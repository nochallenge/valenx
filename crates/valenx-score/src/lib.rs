//! # valenx-score
//!
//! The physics-grounded scoring layer of a biologic-design pipeline: turn a
//! candidate's heterogeneous evidence — interface confidence (ipTM, pLDDT), a
//! docking score, an MM-GBSA-style binding-energy estimate — into **one
//! comparable, ranked value with the components broken out**.
//!
//! ## What
//!
//! - [`energy`] — the textbook molecular-energy primitives, each validated
//!   against its analytic limits: [`energy::coulomb`], [`energy::lennard_jones`]
//!   (and its [`energy::lj_min_distance`]), [`energy::born_solvation`],
//!   [`energy::gb_function`] / [`energy::gb_pair_polar`] (generalized Born), and
//!   [`energy::sasa_nonpolar`].
//! - [`mmgbsa`] — a transparent endpoint binding-energy **decomposition**
//!   ([`mmgbsa::BindingEnergy`]: van der Waals + electrostatic + polar and
//!   nonpolar solvation) and a simple two-bead estimator that assembles one
//!   from the primitives.
//! - [`score`] — [`score::ComparableScore`]: fuse the present evidence channels
//!   ([`score::ScoreComponents`]) into a single `[0, 1]` rank, keeping each
//!   normalized component visible.
//!
//! ## Model
//!
//! The energy primitives are the standard force-field / implicit-solvent terms
//! (Coulomb's law, the 12-6 Lennard-Jones potential, the Born and Still
//! generalized-Born solvation energies, a SASA-linear nonpolar term). The
//! MM-GBSA endpoint approximates a binding free energy as the change in
//! molecular-mechanics energy plus solvation on complex formation. The
//! comparable score maps each evidence channel to `[0, 1]` (confidences pass
//! through; energies and dock scores go through a monotone "more-negative-is-
//! better" transform) and takes a weighted mean over whichever channels are
//! present.
//!
//! ## Honest scope
//!
//! Research/educational grade. The primitives are exact and tested, but the
//! binding-energy estimator is a **simplified single-structure endpoint**, not
//! the real method: a credible MM-GBSA/MM-PBSA ΔG needs an explicit
//! all-atom structure, an MD ensemble, and careful parameterization, and even
//! then correlates only moderately with measured affinity. Treat the ΔG here as
//! an ordering heuristic, and **calibrate** the comparable score against
//! experimental data (e.g. SKEMPI) with [`valenx-calibrate`](https://github.com/nochallenge/valenx)
//! before trusting any probability. Nothing here is a validated affinity or
//! safety prediction.
//!
//! ## Example
//!
//! ```
//! use valenx_score::energy::{lennard_jones, lj_min_distance};
//!
//! // The 12-6 Lennard-Jones well bottoms out at depth -epsilon.
//! let (eps, sigma) = (0.5, 3.4);
//! let r_min = lj_min_distance(sigma);
//! let e = lennard_jones(eps, sigma, r_min).unwrap();
//! assert!((e + eps).abs() < 1e-9);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod energy;
pub mod error;
pub mod mmgbsa;
pub mod score;

pub use error::ScoreError;
pub use mmgbsa::{Bead, BindingEnergy};
pub use score::{ComparableScore, ScoreComponents};
