//! # valenx-marine
//!
//! Box-form ship / boat **hull hydrostatics**: floating displacement and
//! initial transverse stability from a hull's principal dimensions.
//!
//! ## What
//!
//! Given a hull's length, beam, draft, block coefficient and the height of
//! its centre of gravity, this crate computes:
//!
//! - the displaced volume `nabla = Cb * L * B * T` and the displacement
//!   mass `Delta = rho * nabla` (Archimedes), via [`Hull::displaced_volume`]
//!   and [`Hull::displacement_mass`];
//! - the height of the centre of buoyancy above the keel `KB`
//!   ([`Hull::kb`]);
//! - the transverse metacentric radius `BM = I_T / nabla`
//!   ([`Hull::transverse_bm`]); and
//! - the transverse metacentric height `GM = KB + BM - KG`, the standard
//!   initial-stability criterion (`GM > 0` is stable),
//!   ([`Hull::metacentric_height`]).
//!
//! The inverse [`equilibrium_draft`] returns the draft at which a hull of a
//! given mass floats in equilibrium (`Delta` equals the displaced mass).
//!
//! ## Model
//!
//! For a wall-sided hull with a rectangular waterplane `L * B`:
//!
//! - displaced volume `nabla = Cb * L * B * T`, where `Cb` is the block
//!   coefficient (the fraction of the enclosing `L * B * T` box the
//!   immersed hull fills, `0 < Cb <= 1`);
//! - displacement mass `Delta = rho * nabla`, with `rho` the water density
//!   (seawater about 1025 kg/m^3, freshwater about 1000 kg/m^3);
//! - centre of buoyancy above keel `KB = T / 2` (the wall-sided / box
//!   approximation: the centroid of a prismatic immersed volume);
//! - transverse second moment of the waterplane area `I_T = L * B^3 / 12`
//!   (a rectangle about its own centreline);
//! - metacentric radius `BM = I_T / nabla`; and
//! - metacentric height `GM = KB + BM - KG`.
//!
//! Solving `Delta = rho * Cb * L * B * T` for `T` gives the equilibrium
//! draft `T = mass / (rho * Cb * L * B)`.
//!
//! Units are SI and must be self-consistent: metres for `L`, `B`, `T` and
//! `KG`; kg/m^3 for `rho`. Results are m^3 for volume, kg for mass and
//! metres for `KB`, `BM` and `GM`.
//!
//! ## Honest scope
//!
//! Research/educational grade. This is the classic *initial* (small-angle)
//! hydrostatics of a single wall-sided box-form hull on an even keel in
//! calm water. It deliberately omits the real naval-architecture detail a
//! design needs, including but not limited to: the actual hull form and its
//! true `KB` / `I_T` (Morrish or Simpson integration of the offsets rather
//! than the box approximation), the longitudinal metacentre and trim,
//! free-surface corrections, large-angle stability and the righting-arm
//! (GZ) curve, added-weight / free-flooding / damaged stability, and any
//! class / SOLAS / IMO regulatory criteria. It is NOT a substitute for a
//! licensed naval architect or class-approved stability software, and must
//! not be used for design, construction, or any decision affecting safety.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use serde::{Deserialize, Serialize};

/// Standard gravity (m/s^2).
pub const GRAVITY: f64 = 9.80665;
/// Nominal seawater density (kg/m^3).
pub const SEAWATER_DENSITY: f64 = 1025.0;
/// Nominal freshwater density (kg/m^3).
pub const FRESHWATER_DENSITY: f64 = 1000.0;

/// An out-of-domain hydrostatics input.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MarineError {
    /// A quantity that must be finite and strictly positive was not.
    #[error("{quantity} must be finite and positive, got {value}")]
    NonPositive {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },
    /// The block coefficient was outside the physical range `(0, 1]`.
    #[error("block coefficient must be in (0, 1], got {value}")]
    BlockCoefficient {
        /// The offending value.
        value: f64,
    },
    /// A quantity that must be finite (but may be zero) was not.
    #[error("{quantity} must be finite, got {value}")]
    NonFinite {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },
}

/// Return `value` when it is finite and strictly positive, else an error.
fn require_positive(quantity: &'static str, value: f64) -> Result<f64, MarineError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(MarineError::NonPositive { quantity, value })
    }
}

/// Validate a block coefficient lies in the physical range `(0, 1]`.
fn check_block_coefficient(value: f64) -> Result<f64, MarineError> {
    if value.is_finite() && value > 0.0 && value <= 1.0 {
        Ok(value)
    } else {
        Err(MarineError::BlockCoefficient { value })
    }
}

/// A wall-sided box-form hull, validated on construction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Hull {
    /// Waterline length `L` (m).
    pub length_m: f64,
    /// Waterline beam / breadth `B` (m).
    pub beam_m: f64,
    /// Draft `T` — immersed depth (m).
    pub draft_m: f64,
    /// Block coefficient `Cb` in `(0, 1]`.
    pub block_coefficient: f64,
    /// Height of the centre of gravity above the keel `KG` (m).
    pub kg_m: f64,
    /// Water density `rho` (kg/m^3).
    pub water_density: f64,
}

impl Hull {
    /// Build a validated hull. Length, beam, draft and water density must be
    /// finite and positive; the block coefficient must lie in `(0, 1]`; and
    /// `KG` must be finite and non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`MarineError`] when any input is out of its physical domain.
    pub fn new(
        length_m: f64,
        beam_m: f64,
        draft_m: f64,
        block_coefficient: f64,
        kg_m: f64,
        water_density: f64,
    ) -> Result<Self, MarineError> {
        let length_m = require_positive("length", length_m)?;
        let beam_m = require_positive("beam", beam_m)?;
        let draft_m = require_positive("draft", draft_m)?;
        let water_density = require_positive("water density", water_density)?;
        let block_coefficient = check_block_coefficient(block_coefficient)?;
        if !(kg_m.is_finite() && kg_m >= 0.0) {
            return Err(MarineError::NonFinite {
                quantity: "KG",
                value: kg_m,
            });
        }
        Ok(Self {
            length_m,
            beam_m,
            draft_m,
            block_coefficient,
            kg_m,
            water_density,
        })
    }

    /// Displaced volume `nabla = Cb * L * B * T` (m^3).
    pub fn displaced_volume(&self) -> f64 {
        self.block_coefficient * self.length_m * self.beam_m * self.draft_m
    }

    /// Displacement mass `Delta = rho * nabla` (kg).
    pub fn displacement_mass(&self) -> f64 {
        self.water_density * self.displaced_volume()
    }

    /// Displacement in tonnes (1 t = 1000 kg).
    pub fn displacement_tonnes(&self) -> f64 {
        self.displacement_mass() / 1000.0
    }

    /// Buoyant weight `Delta * g` (N).
    pub fn buoyant_weight(&self) -> f64 {
        self.displacement_mass() * GRAVITY
    }

    /// Centre of buoyancy above the keel `KB = T / 2` (m) — the wall-sided
    /// box approximation.
    pub fn kb(&self) -> f64 {
        self.draft_m / 2.0
    }

    /// Transverse second moment of the (rectangular) waterplane area
    /// `I_T = L * B^3 / 12` (m^4).
    pub fn waterplane_inertia_transverse(&self) -> f64 {
        self.length_m * self.beam_m.powi(3) / 12.0
    }

    /// Transverse metacentric radius `BM = I_T / nabla` (m).
    pub fn transverse_bm(&self) -> f64 {
        self.waterplane_inertia_transverse() / self.displaced_volume()
    }

    /// Transverse metacentric height `GM = KB + BM - KG` (m). Positive is
    /// initially (small-angle) stable.
    pub fn metacentric_height(&self) -> f64 {
        self.kb() + self.transverse_bm() - self.kg_m
    }

    /// Whether the hull is initially stable, i.e. `GM > 0`.
    pub fn is_stable(&self) -> bool {
        self.metacentric_height() > 0.0
    }

    /// Compute the full hydrostatics report in one call.
    pub fn hydrostatics(&self) -> Hydrostatics {
        Hydrostatics {
            displaced_volume_m3: self.displaced_volume(),
            displacement_mass_kg: self.displacement_mass(),
            displacement_tonnes: self.displacement_tonnes(),
            kb_m: self.kb(),
            bm_m: self.transverse_bm(),
            gm_m: self.metacentric_height(),
            stable: self.is_stable(),
        }
    }
}

/// A hull's computed initial hydrostatics.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Hydrostatics {
    /// Displaced volume `nabla` (m^3).
    pub displaced_volume_m3: f64,
    /// Displacement mass `Delta` (kg).
    pub displacement_mass_kg: f64,
    /// Displacement (tonnes).
    pub displacement_tonnes: f64,
    /// Centre of buoyancy above the keel `KB` (m).
    pub kb_m: f64,
    /// Transverse metacentric radius `BM` (m).
    pub bm_m: f64,
    /// Transverse metacentric height `GM` (m).
    pub gm_m: f64,
    /// Initial small-angle stability (`GM > 0`).
    pub stable: bool,
}

/// Equilibrium draft (m) at which a wall-sided box-form hull of the given
/// `mass_kg` floats, from `Delta = rho * Cb * L * B * T` solved for `T`:
/// `T = mass / (rho * Cb * L * B)`. The inverse of
/// [`Hull::displacement_mass`].
///
/// # Errors
///
/// Returns [`MarineError`] when any input is out of its physical domain.
pub fn equilibrium_draft(
    length_m: f64,
    beam_m: f64,
    block_coefficient: f64,
    water_density: f64,
    mass_kg: f64,
) -> Result<f64, MarineError> {
    let length_m = require_positive("length", length_m)?;
    let beam_m = require_positive("beam", beam_m)?;
    let water_density = require_positive("water density", water_density)?;
    let mass_kg = require_positive("mass", mass_kg)?;
    let block_coefficient = check_block_coefficient(block_coefficient)?;
    Ok(mass_kg / (water_density * block_coefficient * length_m * beam_m))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= EPS * b.abs().max(1.0)
    }

    #[test]
    fn box_barge_displacement_matches_archimedes() {
        // A pure box (Cb = 1): 120 x 20 x 6 m in seawater.
        let hull = Hull::new(120.0, 20.0, 6.0, 1.0, 8.0, SEAWATER_DENSITY).unwrap();
        assert!(close(hull.displaced_volume(), 14400.0));
        assert!(close(hull.displacement_mass(), 1025.0 * 14400.0));
        assert!(close(hull.displacement_tonnes(), 1025.0 * 14400.0 / 1000.0));
    }

    #[test]
    fn box_barge_metacentre_matches_hand_calc() {
        let hull = Hull::new(120.0, 20.0, 6.0, 1.0, 8.0, SEAWATER_DENSITY).unwrap();
        // KB = T/2 = 3; I_T = L*B^3/12 = 120*8000/12 = 80000; nabla = 14400;
        // BM = 80000/14400 = 5.5556; GM = 3 + 5.5556 - 8 = 0.5556.
        assert!(close(hull.kb(), 3.0));
        assert!(close(hull.waterplane_inertia_transverse(), 80000.0));
        assert!(close(hull.transverse_bm(), 80000.0 / 14400.0));
        assert!(close(
            hull.metacentric_height(),
            3.0 + 80000.0 / 14400.0 - 8.0
        ));
        assert!(hull.is_stable());
    }

    #[test]
    fn raising_kg_past_the_metacentre_makes_it_unstable() {
        // Same hull, but a very high KG (cargo stacked high) -> GM < 0.
        let stable = Hull::new(120.0, 20.0, 6.0, 1.0, 8.0, SEAWATER_DENSITY).unwrap();
        let tippy = Hull::new(120.0, 20.0, 6.0, 1.0, 12.0, SEAWATER_DENSITY).unwrap();
        assert!(stable.is_stable());
        assert!(!tippy.is_stable());
        assert!(tippy.metacentric_height() < 0.0);
    }

    #[test]
    fn equilibrium_draft_inverts_displacement_mass() {
        let hull = Hull::new(80.0, 14.0, 4.5, 0.7, 5.0, SEAWATER_DENSITY).unwrap();
        let mass = hull.displacement_mass();
        let t = equilibrium_draft(
            hull.length_m,
            hull.beam_m,
            hull.block_coefficient,
            hull.water_density,
            mass,
        )
        .unwrap();
        assert!(
            close(t, hull.draft_m),
            "recovered draft {t} != {}",
            hull.draft_m
        );
    }

    #[test]
    fn seawater_floats_higher_than_freshwater() {
        // Same mass: the denser seawater needs less draft than freshwater.
        let mass = 5000000.0;
        let t_sea = equilibrium_draft(60.0, 12.0, 0.8, SEAWATER_DENSITY, mass).unwrap();
        let t_fresh = equilibrium_draft(60.0, 12.0, 0.8, FRESHWATER_DENSITY, mass).unwrap();
        assert!(t_sea < t_fresh);
    }

    #[test]
    fn rejects_out_of_domain_inputs() {
        assert!(Hull::new(0.0, 20.0, 6.0, 1.0, 8.0, SEAWATER_DENSITY).is_err());
        assert!(Hull::new(120.0, 20.0, 6.0, 1.5, 8.0, SEAWATER_DENSITY).is_err());
        assert!(Hull::new(120.0, 20.0, 6.0, 0.0, 8.0, SEAWATER_DENSITY).is_err());
        assert!(Hull::new(120.0, 20.0, 6.0, 1.0, -1.0, SEAWATER_DENSITY).is_err());
        assert!(equilibrium_draft(120.0, 20.0, 1.0, SEAWATER_DENSITY, -1.0).is_err());
    }

    #[test]
    fn hydrostatics_report_is_serde_round_trippable() {
        let hull = Hull::new(100.0, 18.0, 5.5, 0.65, 6.5, SEAWATER_DENSITY).unwrap();
        let h = hull.hydrostatics();
        let json = serde_json::to_string(&h).unwrap();
        let back: Hydrostatics = serde_json::from_str(&json).unwrap();
        // Tolerance compare: serde_json's default float parser is not
        // bit-exact (round-trips to within ~1 ULP), so compare values, not
        // the re-serialized string.
        assert!(close(back.displaced_volume_m3, h.displaced_volume_m3));
        assert!(close(back.displacement_mass_kg, h.displacement_mass_kg));
        assert!(close(back.kb_m, h.kb_m));
        assert!(close(back.bm_m, h.bm_m));
        assert!(close(back.gm_m, h.gm_m));
        assert_eq!(back.stable, h.stable);
    }
}
