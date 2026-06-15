//! Tonicity — comparing two solutions across a water-permeable membrane.
//!
//! **Tonicity** classifies an external solution by how it would change a
//! cell's volume, which is governed by the gradient of **effective**
//! (non-penetrating) osmolytes — the particles the membrane holds back.
//! Water flows toward the side with the higher effective osmolarity.
//!
//! Relative to a reference (e.g. the cell interior), an external
//! solution is:
//!
//! - **Hypotonic**  — lower effective osmolarity than the cell; water
//!   moves *into* the cell (it swells / may lyse).
//! - **Isotonic**   — equal effective osmolarity; **no net water
//!   movement**.
//! - **Hypertonic** — higher effective osmolarity; water moves *out* of
//!   the cell (it shrinks / crenates).
//!
//! ## Tonicity is not osmolarity
//!
//! Osmolarity (see [`crate::vant_hoff::Solution::osmolarity_osmol_per_l`])
//! counts *all* dissolved particles. Tonicity counts only those that
//! **cannot freely cross** the membrane. A solute that equilibrates
//! across the membrane (the classic example is urea, or a cell-permeant
//! cryoprotectant) raises osmolarity but contributes *nothing* to
//! tonicity, because at steady state its concentration is equal on both
//! sides and it exerts no sustained osmotic gradient. This module
//! therefore works on **effective osmolarity** — the osmolarity of the
//! non-penetrating solutes only — and it is the caller's job to pass
//! those figures.

use crate::error::OsmosisError;
use serde::{Deserialize, Serialize};

/// Tonicity of an external solution relative to a reference compartment
/// (conventionally, the inside of a cell).
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Tonicity {
    /// External effective osmolarity **below** the reference: water
    /// enters the cell (it swells).
    Hypotonic,
    /// External effective osmolarity **equal** to the reference: no net
    /// water movement.
    Isotonic,
    /// External effective osmolarity **above** the reference: water
    /// leaves the cell (it shrinks).
    Hypertonic,
}

/// Direction of net water movement implied by a tonicity comparison,
/// stated from the cell's point of view.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum WaterMovement {
    /// Water flows **into** the cell (hypotonic surroundings).
    IntoCell,
    /// **No net** water flow (isotonic surroundings).
    NoNet,
    /// Water flows **out of** the cell (hypertonic surroundings).
    OutOfCell,
}

impl Tonicity {
    /// The net water-movement direction this tonicity produces, from the
    /// cell's frame.
    ///
    /// Water always moves toward higher effective osmolarity, so the
    /// mapping is fixed: hypotonic surroundings draw water *in*,
    /// hypertonic surroundings draw it *out*, isotonic surroundings move
    /// none.
    pub fn water_movement(&self) -> WaterMovement {
        match self {
            Tonicity::Hypotonic => WaterMovement::IntoCell,
            Tonicity::Isotonic => WaterMovement::NoNet,
            Tonicity::Hypertonic => WaterMovement::OutOfCell,
        }
    }
}

/// Classify the tonicity of an `external` solution relative to a
/// `reference` (cell-interior) compartment, comparing **effective
/// osmolarities** (osmol/L of non-penetrating solute).
///
/// `tol` is an absolute osmolarity tolerance (osmol/L): if the two
/// effective osmolarities are within `+/- tol` of each other the result
/// is [`Tonicity::Isotonic`]. `tol` is clamped to be non-negative.
///
/// # Errors
///
/// Returns [`OsmosisError::InvalidParameter`] if either effective
/// osmolarity is non-finite or negative.
pub fn classify_with_tol(
    external_effective_osmolarity: f64,
    reference_effective_osmolarity: f64,
    tol: f64,
) -> Result<Tonicity, OsmosisError> {
    if !external_effective_osmolarity.is_finite() || external_effective_osmolarity < 0.0 {
        return Err(OsmosisError::invalid(
            "external_effective_osmolarity",
            external_effective_osmolarity,
            "must be finite and >= 0",
        ));
    }
    if !reference_effective_osmolarity.is_finite() || reference_effective_osmolarity < 0.0 {
        return Err(OsmosisError::invalid(
            "reference_effective_osmolarity",
            reference_effective_osmolarity,
            "must be finite and >= 0",
        ));
    }
    let tol = tol.max(0.0);
    let diff = external_effective_osmolarity - reference_effective_osmolarity;
    if diff.abs() <= tol {
        Ok(Tonicity::Isotonic)
    } else if diff < 0.0 {
        Ok(Tonicity::Hypotonic)
    } else {
        Ok(Tonicity::Hypertonic)
    }
}

/// [`classify_with_tol`] with a small default tolerance (`1e-9`
/// osmol/L), suitable for exact textbook inputs.
///
/// # Errors
///
/// Propagates [`OsmosisError::InvalidParameter`] from
/// [`classify_with_tol`].
pub fn classify(
    external_effective_osmolarity: f64,
    reference_effective_osmolarity: f64,
) -> Result<Tonicity, OsmosisError> {
    classify_with_tol(
        external_effective_osmolarity,
        reference_effective_osmolarity,
        1e-9,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    // Mammalian cells sit near ~0.300 osmol/L effective osmolarity.
    const CELL: f64 = 0.300;

    #[test]
    fn lower_external_is_hypotonic_water_in() {
        // 0.5% NaCl ~= 0.171 osmol/L < 0.300 -> hypotonic.
        let t = classify(0.171, CELL).unwrap();
        assert_eq!(t, Tonicity::Hypotonic);
        assert_eq!(t.water_movement(), WaterMovement::IntoCell);
    }

    #[test]
    fn higher_external_is_hypertonic_water_out() {
        // 3% NaCl ~= 1.03 osmol/L > 0.300 -> hypertonic.
        let t = classify(1.03, CELL).unwrap();
        assert_eq!(t, Tonicity::Hypertonic);
        assert_eq!(t.water_movement(), WaterMovement::OutOfCell);
    }

    // ---- The required VALIDATE property: isotonic -> no net water --

    #[test]
    fn equal_osmolarity_is_isotonic_no_net_movement() {
        // 0.9% NaCl ~= 0.300 osmol/L matched to the cell -> isotonic.
        let t = classify(CELL, CELL).unwrap();
        assert_eq!(t, Tonicity::Isotonic);
        assert_eq!(t.water_movement(), WaterMovement::NoNet);
    }

    #[test]
    fn isotonic_within_tolerance_band() {
        // A 0.0005 osmol/L mismatch under a 0.001 tol still reads iso.
        let t = classify_with_tol(CELL + 0.0005, CELL, 0.001).unwrap();
        assert_eq!(t, Tonicity::Isotonic);
        assert_eq!(t.water_movement(), WaterMovement::NoNet);
    }

    #[test]
    fn just_outside_tolerance_resolves_to_a_gradient() {
        // +0.01 osmol/L with a tight 1e-9 tol is unambiguously hyper.
        let hyper = classify(CELL + 0.01, CELL).unwrap();
        assert_eq!(hyper, Tonicity::Hypertonic);
        // -0.01 is hypotonic.
        let hypo = classify(CELL - 0.01, CELL).unwrap();
        assert_eq!(hypo, Tonicity::Hypotonic);
    }

    #[test]
    fn water_moves_toward_higher_effective_osmolarity() {
        // Sanity on the direction mapping across all three cases.
        assert_eq!(
            Tonicity::Hypotonic.water_movement(),
            WaterMovement::IntoCell
        );
        assert_eq!(Tonicity::Isotonic.water_movement(), WaterMovement::NoNet);
        assert_eq!(
            Tonicity::Hypertonic.water_movement(),
            WaterMovement::OutOfCell
        );
    }

    #[test]
    fn negative_tolerance_is_clamped_to_zero() {
        // A negative tol must not invert the band; equal values stay iso.
        let t = classify_with_tol(CELL, CELL, -5.0).unwrap();
        assert_eq!(t, Tonicity::Isotonic);
    }

    #[test]
    fn rejects_invalid_osmolarities() {
        assert!(classify(-0.1, CELL).is_err());
        assert!(classify(CELL, -0.1).is_err());
        assert!(classify(f64::NAN, CELL).is_err());
        assert!(classify(CELL, f64::INFINITY).is_err());
    }
}
