//! Docking analysis, validation, library filtering, I/O and covalent
//! docking.
//!
//! This module collects the post-docking and screening-support
//! features that are not search or scoring proper:
//!
//! - [`pocket`] — geometric binding-pocket detection (an fpocket-class
//!   grid / sphere method).
//! - [`fingerprint`] — a protein-ligand interaction fingerprint
//!   (H-bond, hydrophobic, π-stacking, salt-bridge, halogen contacts).
//! - [`validate`] — pose-RMSD computation and redocking validation
//!   (the success-rate metric used to benchmark a docking protocol).
//! - [`redock_bench`] — inline canonical PDB co-crystal cases (1HVR,
//!   3PTB, 1STP) and a [`redock_bench::run_canonical_benchmark`]
//!   one-call validation pass over them.
//! - [`pharmacophore`] — pharmacophore-based library screening,
//!   reusing [`valenx_cheminf`]'s pharmacophore features.
//! - [`admet`] — ADMET-lite library filtering for screening, reusing
//!   [`valenx_cheminf`]'s descriptors, Lipinski / Veber rules and the
//!   structural-alert (PAINS-lite) count.
//! - [`io`] — PDBQT-style docking-result reading / writing and a
//!   results table.
//! - [`covalent`] — covalent-docking driver, anchoring the ligand at a
//!   reactive receptor atom.

pub mod admet;
pub mod covalent;
pub mod fingerprint;
pub mod io;
pub mod pharmacophore;
pub mod pocket;
pub mod redock_bench;
pub mod validate;

pub use admet::{admet_filter, AdmetFilter, AdmetVerdict};
pub use covalent::{covalent_dock, CovalentParams, CovalentResult};
pub use fingerprint::{interaction_fingerprint, InteractionFingerprint, InteractionKind};
pub use io::{read_docking_result, write_docking_result, DockingResultTable};
pub use pharmacophore::{pharmacophore_screen, PharmacophoreQuery};
pub use pocket::{detect_pockets, Pocket};
pub use redock_bench::{
    canonical_cases, hiv1_protease_1hvr, run_canonical_benchmark, streptavidin_1stp, trypsin_3ptb,
};
pub use validate::{redock_success_rate, RedockOutcome};
