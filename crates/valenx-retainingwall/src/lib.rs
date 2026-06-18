//! # valenx-retainingwall
//!
//! Lateral earth pressure on a retaining wall by **Rankine's theory**.
//!
//! ## What
//!
//! Given a soil's internal friction angle `phi` and unit weight `gamma`,
//! this crate computes the classic Rankine answers for a smooth vertical
//! wall retaining dry, cohesionless, level backfill:
//!
//! - the **active** and **passive** earth-pressure coefficients
//!   `Ka = tan^2(45 deg - phi/2)` and `Kp = tan^2(45 deg + phi/2) = 1/Ka`;
//! - the **lateral pressure** at any depth, `sigma = K * gamma * z`,
//!   which grows linearly with depth;
//! - the **resultant thrust** per unit length of wall,
//!   `P = 1/2 * K * gamma * H^2`, together with its **line of action** at
//!   `H/3` above the base;
//! - the **inverse** of that thrust — the wall height
//!   `H = sqrt(2 P / (K gamma))` that produces a target resultant
//!   ([`SoilProfile::height_for_active_thrust`] /
//!   [`SoilProfile::height_for_passive_thrust`]).
//!
//! ```
//! use valenx_retainingwall::SoilProfile;
//!
//! // Medium-dense sand: phi = 30 deg, gamma = 18 kN/m^3.
//! let soil = SoilProfile::new(30.0, 18.0).expect("valid soil");
//!
//! // Active coefficient Ka = tan^2(30 deg) = 1/3.
//! assert!((soil.ka() - 1.0 / 3.0).abs() < 1e-9);
//!
//! // Thrust on a 5 m wall: Pa = 1/2 * (1/3) * 18 * 25 = 75 kN/m at H/3.
//! let thrust = soil.active_thrust(5.0).expect("valid height");
//! assert!((thrust.resultant - 75.0).abs() < 1e-9);
//! assert!((thrust.line_of_action - 5.0 / 3.0).abs() < 1e-9);
//! ```
//!
//! ## Model
//!
//! Rankine's theory treats the soil behind a frictionless vertical wall
//! as a cohesionless granular mass on the verge of plastic failure. The
//! horizontal-to-vertical effective-stress ratio at failure is the
//! earth-pressure coefficient `K`; for a level surface it depends only on
//! `phi`. Active conditions (wall moving away from the soil) mobilise the
//! minimum ratio `Ka < 1`; passive conditions (wall pushed into the soil)
//! mobilise the maximum ratio `Kp = 1/Ka > 1`. The vertical effective
//! stress is hydrostatic in the soil, `sigma_v = gamma * z`, so the
//! lateral pressure profile is the triangle `K * gamma * z`, and its
//! resultant is the triangle's area `1/2 * K * gamma * H^2` acting at the
//! centroid `H/3` above the base. Inverting the resultant gives the wall
//! height for a target thrust, `H = sqrt(2 P / (K gamma))`. See
//! [`rankine`] for the full derivation and the [`rankine::SoilProfile`]
//! API.
//!
//! ## Honest scope
//!
//! Research/educational grade. This crate implements only the
//! introductory textbook closed form: a smooth (frictionless) vertical
//! wall, dry cohesionless soil (`c = 0`), and a horizontal backfill with
//! no surcharge or water table. It deliberately does **not** model
//! cohesion, wall friction (the Coulomb wedge), sloping or surcharged
//! backfill, pore-water pressure, layered or partially saturated soils,
//! seismic (Mononobe-Okabe) earth pressure, or any wall-stability limit
//! state (sliding, overturning, bearing capacity, global stability). It
//! is **not** a clinical/medical tool and **not** a production
//! geotechnical-engineering or structural-design tool; do not use it for
//! the design of real retaining structures.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod rankine;

pub use error::{ErrorCategory, RetainingWallError};
pub use rankine::{SoilProfile, Thrust};
