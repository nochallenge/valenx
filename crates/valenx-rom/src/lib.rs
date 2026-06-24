//! # valenx-rom — reduced-order & surrogate modelling
//!
//! Turn time-resolved state data into fast, low-dimensional models. Given a
//! [`Snapshots`] matrix whose **columns are the system state at successive
//! times**, this crate provides the four workhorse reduced-order modelling
//! (ROM) techniques, all built on the workspace [`nalgebra`] linear-algebra
//! stack (SVD / least-squares / complex eigen-solve) and
//! [`num_complex`] for spectra:
//!
//! - **POD** — [`PodBasis`]: proper orthogonal decomposition. The energy-optimal
//!   orthonormal basis via the SVD of the snapshot matrix, truncated by
//!   cumulative singular-value energy. [`PodBasis::project`] /
//!   [`PodBasis::reconstruct`] / [`PodBasis::reconstruction_error`].
//! - **DMD** — [`Dmd`]: dynamic mode decomposition (standard & exact). The
//!   best-fit linear operator's spectrum — complex [`Dmd::eigenvalues`],
//!   spatial [`Dmd::modes`], and continuous-time [`Dmd::growth_rates`] /
//!   [`Dmd::frequencies`] — with optional POD rank truncation.
//! - **OpInf** — [`OpInfModel`]: operator inference. Non-intrusive least-squares
//!   identification of reduced linear (and optional quadratic) operators from
//!   projected snapshots, with reduced-state RK4 time-stepping
//!   ([`OpInfModel::step`] / [`OpInfModel::rollout`]).
//! - **POD–Galerkin** — [`galerkin::galerkin_operator`]: intrusive projection of
//!   a known full-order linear operator onto a POD basis.
//!
//! ## Fail-loud contract
//!
//! Valenx computes numbers people make decisions from, so every routine here
//! returns a [`RomError`] rather than a plausible-but-wrong value on bad input:
//! empty data, a numerically rank-deficient snapshot set, a dimension mismatch,
//! a non-finite entry, an out-of-range energy tolerance, an impossible
//! truncation rank, or too few time samples. Nothing panics on user input.
//!
//! ## Honest scope
//!
//! These are textbook, deterministic methods. POD/DMD are exact linear-algebra
//! constructions (validated to machine precision against analytic fields and a
//! known complex eigenvalue in the test suite). OpInf is a *data-driven* fit —
//! its accuracy depends on snapshot quality, the finite-difference derivative
//! estimate (O(dt²) here), the regularisation (see [`opinf::OpInfConfig::ridge`]),
//! and how well the chosen operator structure matches the true dynamics; it can
//! be ill-conditioned and is not a turnkey black box. DMD rank truncation
//! smaller than the data rank yields *approximate* eigenvalues, not exact ones.
//! Calibrate energy tolerances, ranks, and regularisation against your problem.
//!
//! ## Example — POD round trip
//!
//! ```
//! use valenx_rom::{Snapshots, PodBasis};
//! use nalgebra::DMatrix;
//!
//! // Two columns spanning a 2-D plane in 3-D state space.
//! let x = DMatrix::from_row_slice(3, 4, &[
//!     1.0, 2.0, 3.0, 4.0,
//!     0.5, 1.0, 1.5, 2.0,
//!     0.0, 0.0, 0.0, 0.0,
//! ]);
//! let snaps = Snapshots::from_matrix(x).unwrap();
//! let basis = PodBasis::fit(&snaps, 0.999).unwrap();
//! assert_eq!(basis.rank(), 1);                 // data is rank 1
//! assert!(basis.reconstruction_error(&snaps).unwrap() < 1e-10);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod dmd;
pub mod error;
pub mod galerkin;
pub mod opinf;
pub mod pod;
pub mod snapshots;

pub use dmd::{Dmd, DmdVariant};
pub use error::RomError;
pub use galerkin::{galerkin_affine, galerkin_operator};
pub use opinf::{forward_difference, kron_unique, OpInfConfig, OpInfModel};
pub use pod::PodBasis;
pub use snapshots::Snapshots;
