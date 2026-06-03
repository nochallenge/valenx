//! Cryo-EM reconstruction — classical signal processing (non-ML).
//!
//! **Roadmap features 23-29.** Single-particle cryo-electron
//! microscopy reconstructs a 3-D density map of a molecule from tens
//! of thousands of noisy 2-D projection images of individual copies
//! of it, frozen in random orientations. The classical reconstruction
//! pipeline — the signal-processing core of RELION, EMAN2 and cisTEM
//! — is a sequence of well-defined algorithms, *none* of which is a
//! neural network:
//!
//! - **MRC I/O + data model** ([`mrc`]) — read / write the MRC
//!   density-map format; a particle-stack / micrograph model.
//! - **CTF** ([`ctf`]) — the contrast-transfer function model and its
//!   estimation from a power spectrum.
//! - **Particle picking** ([`picking`]) — find particles in a
//!   micrograph by template / blob correlation.
//! - **2D class averaging** ([`classify`]) — align particles in 2-D
//!   and average them into high-SNR class averages.
//! - **3D reconstruction** ([`reconstruct`]) — build a 3-D map from
//!   oriented 2-D projections by weighted back-projection / direct
//!   Fourier inversion.
//! - **Projection-matching refinement** ([`refine_em`]) — iteratively
//!   assign particle orientations by matching reprojections of the
//!   current map, then reconstruct, then repeat.
//! - **Resolution estimation** ([`fsc`]) — Fourier shell correlation
//!   and the gold-standard half-map criterion.
//!
//! **Honest note.** This is the *classical* reconstruction core. It
//! is real and correct — back-projection genuinely inverts the
//! projection geometry, FSC genuinely measures resolution. It is *not*
//! the full RELION pipeline: RELION wraps these steps in a regularised
//! maximum-likelihood (empirical-Bayes) framework, runs on a GPU
//! cluster, and the modern particle pickers are neural networks (those
//! stay adapter-only). This module is the textbook signal-processing
//! backbone, honestly scoped.

pub mod classify;
pub mod ctf;
pub mod fsc;
pub mod mrc;
pub mod picking;
pub mod reconstruct;
pub mod refine_em;

pub use classify::{class_averages, ClassAverageResult};
pub use ctf::{estimate_ctf, Ctf, CtfEstimate};
pub use fsc::{fourier_shell_correlation, gold_standard_resolution, FscCurve};
pub use mrc::{Image2d, ParticleStack, Volume3d};
pub use picking::{pick_particles, ParticlePick};
pub use reconstruct::{reconstruct_3d, Projection, ReconstructionResult};
pub use refine_em::{projection_matching, RefinementResult};
