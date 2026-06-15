//! # valenx-insulation
//!
//! Steady-state building-envelope heat-loss models: thermal
//! resistance, surface films, composite walls, U-values, and
//! conductive heat-loss rate.
//!
//! ## What
//!
//! Given the materials that make up a wall, this crate answers three
//! coupled questions a building-physics textbook asks:
//!
//! 1. How much does each element resist heat flow? — the area-specific
//!    resistance ("R-value") of a solid [`Layer`] or a convective
//!    [`SurfaceFilm`].
//! 2. How well does the whole wall keep heat in? — the total series
//!    resistance and the U-value of a [`CompositeWall`].
//! 3. How fast does heat leak out? — the steady-state heat-loss rate
//!    `Q` for a given wall area and inside/outside temperature
//!    difference.
//!
//! ## Model
//!
//! Heat flow is treated as one-dimensional and steady. Under that
//! assumption the thermal resistances of a wall's elements add in
//! series, and the governing relations are the standard closed forms:
//!
//! - Conduction through a solid layer: `R = L / k`, with thickness `L`
//!   (m) and conductivity `k` (`W/(m.K)`); units `m^2.K/W`.
//! - Convective/radiative surface film: `R_film = 1 / h`, with film
//!   coefficient `h` (`W/(m^2.K)`).
//! - Composite wall (series): `R_total = R_si + sum(L_i / k_i) + R_se`.
//! - Thermal transmittance: `U = 1 / R_total` (`W/(m^2.K)`).
//! - Heat-loss rate: `Q = U * A * dT` (W), with area `A` (m^2) and
//!   temperature difference `dT` (K).
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form,
//! one-dimensional steady-state models. They deliberately ignore
//! two-dimensional thermal bridging, transient/dynamic storage,
//! moisture and vapour transport, air infiltration, solar gains, and
//! the temperature dependence of conductivity. This crate is NOT a
//! clinical/medical/production engineering tool: do not use it for
//! building-code compliance energy modelling, HVAC equipment sizing,
//! or any safety-critical decision. For those, use a validated,
//! standards-conformant tool and a qualified engineer.
//!
//! ## Example
//!
//! ```
//! use valenx_insulation::{CompositeWall, Layer, SurfaceFilm};
//!
//! // A simple insulated wall between standard interior/exterior films.
//! let wall = CompositeWall::builder()
//!     .interior_film(SurfaceFilm::interior_default())
//!     .layer(Layer::new(0.16, 0.04).unwrap()) // 160 mm of k = 0.04 insulation
//!     .exterior_film(SurfaceFilm::exterior_default())
//!     .build()
//!     .unwrap();
//!
//! // R = 0.13 + 4.0 + 0.04 = 4.17 m^2.K/W.
//! assert!((wall.total_resistance() - 4.17).abs() < 1e-3);
//!
//! // Heat loss through 10 m^2 at a 20 K difference.
//! let q = wall.heat_loss(10.0, 20.0).unwrap();
//! assert!(q > 47.0 && q < 49.0);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod thermal;
pub mod wall;

pub use error::{ErrorCategory, InsulationError};
pub use thermal::{Layer, SurfaceFilm};
pub use wall::{CompositeWall, CompositeWallBuilder};
