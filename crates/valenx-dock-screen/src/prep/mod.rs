//! Receptor and ligand preparation.
//!
//! Before a docking search can run, three things have to be settled:
//!
//! - **Where** to dock — the [`gridbox`] module defines the cubic /
//!   rectangular search volume (centre + size + spacing).
//! - **What can move** — the [`torsion`] module extracts a ligand's
//!   torsion tree (a rigid root plus a set of rotatable branches), and
//!   [`flex`] selects which receptor side chains are allowed to flex.
//! - **How atoms interact** — the [`protonate`] module assigns a
//!   protonation state at a target pH and Gasteiger-style partial
//!   charges, reusing [`valenx_cheminf`]'s PEOE solver.

pub mod flex;
pub mod gridbox;
pub mod protonate;
pub mod torsion;

pub use flex::{FlexibleSidechain, FlexSelection};
pub use gridbox::GridBox;
pub use protonate::{ChargeModel, ProtonationResult};
pub use torsion::{RotatableBond, TorsionTree};
