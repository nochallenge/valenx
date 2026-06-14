//! # valenx-buckling
//!
//! Euler elastic column-buckling calculator: given a slender column's
//! material, cross-section, length, and end restraints, compute the
//! load at which it suddenly bows sideways and the stress and
//! slenderness that go with it.
//!
//! ## What
//!
//! Build a validated [`Column`] from Young's modulus `E`, the smallest
//! second moment of area `I`, the unsupported length `L`, the
//! cross-sectional area `A`, and an [`EndCondition`], then read back:
//!
//! - [`Column::critical_load`] — the **Euler critical load** `P_cr`,
//! - [`Column::critical_stress`] — the **critical stress** `P_cr / A`,
//! - [`Column::slenderness_ratio`] — the **slenderness** `K L / r`,
//! - plus [`Column::effective_length`] (`K L`) and
//!   [`Column::radius_of_gyration`] (`r = sqrt(I / A)`).
//!
//! ```
//! use valenx_buckling::{Column, EndCondition};
//!
//! // A 3 m round steel strut, pinned at both ends.
//! let col = Column::new(
//!     200.0e9,   // E  [Pa]
//!     4.909e-6,  // I  [m^4]   (50 mm dia solid round)
//!     3.0,       // L  [m]
//!     3.142e-3,  // A  [m^2]
//!     EndCondition::PinnedPinned,
//! )
//! .expect("valid column");
//!
//! // P_cr = pi^2 E I / (K L)^2  ~=  1.08 MN.
//! let pcr = col.critical_load();
//! assert!((pcr - 1.0768e6).abs() < 2.0e3);
//!
//! // Clamping both ends (fixed-fixed, K = 0.5) carries 4x as much.
//! let clamped = Column::new(200.0e9, 4.909e-6, 3.0, 3.142e-3, EndCondition::FixedFixed)
//!     .expect("valid column");
//! assert!((clamped.critical_load() / pcr - 4.0).abs() < 1e-9);
//! ```
//!
//! ## Model
//!
//! The single governing relation is the **Euler buckling formula** for
//! a perfectly straight, axially loaded, linear-elastic prismatic
//! column:
//!
//! ```text
//! P_cr = pi^2 E I / (K L)^2
//! ```
//!
//! The **effective-length factor** `K` encodes the four classical end
//! conditions (theoretical AISC / Hibbeler values):
//!
//! | End condition | `K`  |
//! |---------------|------|
//! | Pinned-pinned | 1.0  |
//! | Fixed-free    | 2.0  |
//! | Fixed-fixed   | 0.5  |
//! | Fixed-pinned  | 0.7  |
//!
//! `K L` is the **effective length** (the equivalent pinned-pinned
//! span), so `P_cr` scales as `1 / K^2`: a fixed-fixed column
//! (`K = 0.5`) carries `4x`, and a fixed-free column (`K = 2`) only
//! `1/4x`, the load of an otherwise identical pinned-pinned column. The
//! derived quantities are the **critical stress** `sigma_cr = P_cr / A`
//! and the **slenderness ratio** `K L / r` with radius of gyration
//! `r = sqrt(I / A)`; equivalently `sigma_cr = pi^2 E / (K L / r)^2`, so
//! the Euler stress falls off as `1 / slenderness^2`.
//!
//! ## Honest scope
//!
//! Research/educational grade. This implements the **textbook
//! closed-form elastic (Euler) buckling model** and nothing more — it
//! is **NOT a clinical/medical tool and NOT a production
//! structural-engineering design tool.** In particular it deliberately
//! does *not* model:
//!
//! - **Inelastic / intermediate columns** — no yield cap and no
//!   Johnson (parabolic) or tangent-modulus transition, so for a stocky
//!   column the reported elastic `sigma_cr` can exceed the material
//!   yield strength, where a real column would crush or yield first.
//! - **Imperfections** — initial crookedness, load eccentricity, and
//!   residual stresses (the Perry-Robertson / real-column knockdowns)
//!   are ignored; the column is treated as perfectly straight and
//!   concentrically loaded.
//! - **Code safety factors** — no AISC/Eurocode design factors,
//!   load-and-resistance factoring, or local/torsional/lateral-torsional
//!   buckling of thin-walled sections; only flexural (Euler) buckling
//!   of a prismatic member about a single principal axis.
//!
//! Use it to learn and to sanity-check, not to certify a structure.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod column;
pub mod error;

pub use column::{Column, EndCondition};
pub use error::BucklingError;
