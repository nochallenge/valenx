//! Allowable working stresses for the rivet and plate material.
//!
//! A riveted joint mixes three independent strength checks, each driven
//! by a different allowable stress:
//!
//! - the **rivet** shank resists *shear* across its cross-section,
//! - the **plate** (and rivet) resist *bearing* (crushing) on the
//!   projected contact area, and
//! - the **plate** resists *tension* across its net section.
//!
//! [`Allowables`] groups the three permissible stresses so a joint
//! evaluation has a single source for "how hard may each part work".
//! All stresses are in pascals (`N/m²`).

use crate::error::{Result, RivetError};
use serde::{Deserialize, Serialize};

/// Permissible working stresses for a riveted joint, in pascals.
///
/// Each field is the stress at which the corresponding failure mode is
/// taken to occur — for a *working-stress* (allowable) design these are
/// the code allowables; for an *ultimate* check they are the ultimate
/// strengths. The calculator treats them purely as the stress that
/// multiplies the relevant area to give a failure load, so the same
/// type serves either interpretation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Allowables {
    /// Permissible **shear** stress in the rivet shank, `τ` (Pa).
    pub shear: f64,
    /// Permissible **bearing** (crushing) stress on the plate / rivet
    /// projected contact area, `σ_b` (Pa).
    pub bearing: f64,
    /// Permissible **tensile** stress in the plate, `σ_t` (Pa).
    pub tension: f64,
}

impl Allowables {
    /// Build a validated set of allowable stresses.
    ///
    /// # Errors
    ///
    /// Returns [`RivetError::NotPositive`] if any of `shear`, `bearing`
    /// or `tension` is not finite and strictly positive.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_rivet::Allowables;
    ///
    /// // A common mild-steel rivet/plate working-stress set (MPa):
    /// // τ = 80, σ_b = 160, σ_t = 100.
    /// let a = Allowables::new(80.0e6, 160.0e6, 100.0e6).unwrap();
    /// assert!((a.shear - 80.0e6).abs() < 1.0);
    /// ```
    pub fn new(shear: f64, bearing: f64, tension: f64) -> Result<Self> {
        Ok(Self {
            shear: RivetError::require_positive("shear", shear)?,
            bearing: RivetError::require_positive("bearing", bearing)?,
            tension: RivetError::require_positive("tension", tension)?,
        })
    }
}
