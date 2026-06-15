//! Starling's equation — transcapillary fluid exchange.
//!
//! The volumetric flux of fluid across a capillary wall (the
//! semipermeable membrane separating plasma from the interstitium) is
//!
//! ```text
//! Jv = Kf * ( (Pc - Pi) - sigma * (pi_c - pi_i) )
//! ```
//!
//! with
//!
//! - `Pc`  — capillary (plasma) **hydrostatic** pressure,
//! - `Pi`  — interstitial **hydrostatic** pressure,
//! - `pi_c` — capillary (plasma) **oncotic** (colloid osmotic) pressure,
//! - `pi_i` — interstitial **oncotic** pressure,
//! - `sigma` — Staverman **reflection coefficient** (`0`..=`1`; `1` =
//!   wall perfectly impermeable to the colloid, `0` = freely permeable),
//! - `Kf`  — **filtration coefficient** (membrane hydraulic
//!   conductivity times surface area; `>= 0`).
//!
//! The term `(Pc - Pi)` is the net **hydrostatic** pressure pushing
//! fluid *out* of the capillary; `sigma * (pi_c - pi_i)` is the net
//! **oncotic** pressure pulling it *back in*. Their difference is the
//! **net driving pressure** (NDP); multiplying by `Kf` gives the flux.
//!
//! Sign convention: **positive `Jv` is filtration** (fluid leaving the
//! capillary into the tissue); **negative `Jv` is reabsorption** (fluid
//! returning to the capillary).
//!
//! All pressures share one unit (this crate does not assume which —
//! mmHg is the physiological convention). `Kf` then carries
//! `volume / (time * pressure)`, and `Jv` comes out in `volume / time`.

use crate::error::OsmosisError;
use serde::{Deserialize, Serialize};

/// The six Starling parameters for one capillary segment.
///
/// Build via [`StarlingParams::new`], which validates `sigma` and `Kf`.
/// The four pressures may be any finite value (interstitial hydrostatic
/// pressure, for instance, is slightly **sub-atmospheric** — i.e.
/// negative on a gauge scale — in many tissues).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StarlingParams {
    /// Capillary hydrostatic pressure `Pc`.
    pub capillary_hydrostatic: f64,
    /// Interstitial hydrostatic pressure `Pi`.
    pub interstitial_hydrostatic: f64,
    /// Capillary (plasma) oncotic pressure `pi_c`.
    pub capillary_oncotic: f64,
    /// Interstitial oncotic pressure `pi_i`.
    pub interstitial_oncotic: f64,
    /// Staverman reflection coefficient `sigma`, in `0..=1`.
    pub reflection_sigma: f64,
    /// Filtration coefficient `Kf`, `>= 0`.
    pub filtration_kf: f64,
}

/// Which way fluid is moving across the membrane.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum FluxDirection {
    /// `Jv > 0` — fluid leaves the capillary for the interstitium.
    Filtration,
    /// `Jv < 0` — fluid returns from the interstitium to the capillary.
    Reabsorption,
    /// `Jv == 0` — the hydrostatic and oncotic gradients balance.
    Equilibrium,
}

impl StarlingParams {
    /// Build a validated [`StarlingParams`].
    ///
    /// # Errors
    ///
    /// Returns [`OsmosisError::InvalidParameter`] if any pressure is
    /// non-finite, if `reflection_sigma` is outside `0..=1` (or
    /// non-finite), or if `filtration_kf` is negative (or non-finite).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        capillary_hydrostatic: f64,
        interstitial_hydrostatic: f64,
        capillary_oncotic: f64,
        interstitial_oncotic: f64,
        reflection_sigma: f64,
        filtration_kf: f64,
    ) -> Result<Self, OsmosisError> {
        for (name, v) in [
            ("capillary_hydrostatic", capillary_hydrostatic),
            ("interstitial_hydrostatic", interstitial_hydrostatic),
            ("capillary_oncotic", capillary_oncotic),
            ("interstitial_oncotic", interstitial_oncotic),
        ] {
            if !v.is_finite() {
                return Err(OsmosisError::invalid(name, v, "must be finite"));
            }
        }
        if !reflection_sigma.is_finite() || !(0.0..=1.0).contains(&reflection_sigma) {
            return Err(OsmosisError::invalid(
                "reflection_sigma",
                reflection_sigma,
                "must be finite and in [0, 1]",
            ));
        }
        if !filtration_kf.is_finite() || filtration_kf < 0.0 {
            return Err(OsmosisError::invalid(
                "filtration_kf",
                filtration_kf,
                "must be finite and >= 0",
            ));
        }
        Ok(Self {
            capillary_hydrostatic,
            interstitial_hydrostatic,
            capillary_oncotic,
            interstitial_oncotic,
            reflection_sigma,
            filtration_kf,
        })
    }

    /// Net **hydrostatic** pressure gradient `Pc - Pi` (outward-positive).
    pub fn net_hydrostatic(&self) -> f64 {
        self.capillary_hydrostatic - self.interstitial_hydrostatic
    }

    /// Net **oncotic** pressure gradient `pi_c - pi_i` (the inward pull,
    /// before applying `sigma`).
    pub fn net_oncotic(&self) -> f64 {
        self.capillary_oncotic - self.interstitial_oncotic
    }

    /// Net driving pressure
    /// `NDP = (Pc - Pi) - sigma * (pi_c - pi_i)`.
    ///
    /// Positive favours filtration, negative favours reabsorption. This
    /// is `Jv / Kf`, so its sign alone (independent of `Kf`) decides the
    /// flux direction.
    pub fn net_driving_pressure(&self) -> f64 {
        self.net_hydrostatic() - self.reflection_sigma * self.net_oncotic()
    }

    /// Net transcapillary flux
    /// `Jv = Kf * NDP`.
    ///
    /// Positive = filtration (out of the capillary); negative =
    /// reabsorption (into the capillary).
    pub fn net_filtration(&self) -> f64 {
        self.filtration_kf * self.net_driving_pressure()
    }

    /// Classify the flux direction from the sign of the net driving
    /// pressure.
    ///
    /// Uses an absolute tolerance `tol` (in the shared pressure unit) so
    /// a driving pressure within `+/- tol` of zero reads as
    /// [`FluxDirection::Equilibrium`]. `tol` is clamped to be
    /// non-negative.
    pub fn direction_with_tol(&self, tol: f64) -> FluxDirection {
        let tol = tol.max(0.0);
        let ndp = self.net_driving_pressure();
        if ndp > tol {
            FluxDirection::Filtration
        } else if ndp < -tol {
            FluxDirection::Reabsorption
        } else {
            FluxDirection::Equilibrium
        }
    }

    /// [`direction_with_tol`](Self::direction_with_tol) with a small
    /// default tolerance (`1e-9`) suitable for exact textbook inputs.
    pub fn direction(&self) -> FluxDirection {
        self.direction_with_tol(1e-9)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // Classic arteriolar-end capillary numbers (mmHg), Guyton-style:
    //   Pc = 35, Pi = -2 (slightly sub-atmospheric),
    //   pi_c = 28, pi_i = 5, sigma = 1, Kf = 1 (unit).
    // NDP = (35 - (-2)) - 1*(28 - 5) = 37 - 23 = +14 mmHg -> filtration.
    fn arterial_end() -> StarlingParams {
        StarlingParams::new(35.0, -2.0, 28.0, 5.0, 1.0, 1.0).unwrap()
    }

    // Venular end: capillary hydrostatic pressure has dropped.
    //   Pc = 15, rest unchanged.
    // NDP = (15 - (-2)) - (28 - 5) = 17 - 23 = -6 mmHg -> reabsorption.
    fn venular_end() -> StarlingParams {
        StarlingParams::new(15.0, -2.0, 28.0, 5.0, 1.0, 1.0).unwrap()
    }

    #[test]
    fn arterial_end_filters_outward() {
        let p = arterial_end();
        assert!((p.net_driving_pressure() - 14.0).abs() < EPS);
        assert!(p.net_filtration() > 0.0);
        assert_eq!(p.direction(), FluxDirection::Filtration);
    }

    #[test]
    fn venular_end_reabsorbs() {
        let p = venular_end();
        assert!((p.net_driving_pressure() - (-6.0)).abs() < EPS);
        assert!(p.net_filtration() < 0.0);
        assert_eq!(p.direction(), FluxDirection::Reabsorption);
    }

    // ---- The required VALIDATE property: sign flips when oncotic
    //      pull exceeds the hydrostatic push -------------------------

    #[test]
    fn sign_flips_when_oncotic_exceeds_hydrostatic() {
        // Hold the oncotic gradient fixed (sigma*(pi_c - pi_i) = 23) and
        // sweep the hydrostatic push across it.
        let oncotic_pull = 23.0;

        // Hydrostatic push 30 > 23 -> net OUTWARD (filtration, +).
        let push_big = StarlingParams::new(30.0, 0.0, 28.0, 5.0, 1.0, 1.0).unwrap();
        assert!(push_big.net_filtration() > 0.0);
        assert_eq!(push_big.direction(), FluxDirection::Filtration);

        // Hydrostatic push 10 < 23 -> net INWARD (reabsorption, -).
        let push_small = StarlingParams::new(10.0, 0.0, 28.0, 5.0, 1.0, 1.0).unwrap();
        assert!(push_small.net_filtration() < 0.0);
        assert_eq!(push_small.direction(), FluxDirection::Reabsorption);

        // The crossover sits exactly where hydrostatic == oncotic pull.
        let balanced = StarlingParams::new(oncotic_pull, 0.0, 28.0, 5.0, 1.0, 1.0).unwrap();
        assert!(balanced.net_driving_pressure().abs() < EPS);
        assert_eq!(balanced.direction(), FluxDirection::Equilibrium);
    }

    #[test]
    fn equilibrium_when_gradients_balance() {
        // Pc - Pi = pi_c - pi_i, sigma = 1 -> NDP = 0.
        let p = StarlingParams::new(28.0, 0.0, 33.0, 5.0, 1.0, 5.0).unwrap();
        assert!(p.net_driving_pressure().abs() < EPS);
        assert!(p.net_filtration().abs() < EPS);
        assert_eq!(p.direction(), FluxDirection::Equilibrium);
    }

    #[test]
    fn flux_scales_linearly_with_kf() {
        let ndp_only = StarlingParams::new(35.0, -2.0, 28.0, 5.0, 1.0, 1.0).unwrap();
        let kf_ten = StarlingParams::new(35.0, -2.0, 28.0, 5.0, 1.0, 10.0).unwrap();
        // Same NDP, 10x Kf -> exactly 10x flux.
        assert!((kf_ten.net_filtration() - 10.0 * ndp_only.net_filtration()).abs() < EPS);
    }

    #[test]
    fn zero_kf_gives_zero_flux_regardless_of_gradient() {
        let p = StarlingParams::new(35.0, -2.0, 28.0, 5.0, 1.0, 0.0).unwrap();
        assert!(p.net_driving_pressure().abs() > 1.0); // gradient present
        assert!(p.net_filtration().abs() < EPS); // but no conductance
    }

    #[test]
    fn reflection_coefficient_scales_the_oncotic_term() {
        // sigma = 0 -> the wall ignores the colloid; only hydrostatics
        // drive flow, so NDP = Pc - Pi regardless of oncotic gradient.
        let leaky = StarlingParams::new(20.0, 0.0, 28.0, 5.0, 0.0, 1.0).unwrap();
        assert!((leaky.net_driving_pressure() - 20.0).abs() < EPS);
        assert_eq!(leaky.direction(), FluxDirection::Filtration);

        // sigma = 0.5 -> half the oncotic pull (23 -> 11.5):
        // NDP = 20 - 11.5 = 8.5.
        let half = StarlingParams::new(20.0, 0.0, 28.0, 5.0, 0.5, 1.0).unwrap();
        assert!((half.net_driving_pressure() - 8.5).abs() < EPS);

        // sigma = 1 -> full oncotic pull: NDP = 20 - 23 = -3 -> reabsorb.
        let tight = StarlingParams::new(20.0, 0.0, 28.0, 5.0, 1.0, 1.0).unwrap();
        assert!((tight.net_driving_pressure() - (-3.0)).abs() < EPS);
        assert_eq!(tight.direction(), FluxDirection::Reabsorption);
    }

    #[test]
    fn validation_rejects_out_of_range() {
        assert!(StarlingParams::new(35.0, -2.0, 28.0, 5.0, 1.5, 1.0).is_err()); // sigma>1
        assert!(StarlingParams::new(35.0, -2.0, 28.0, 5.0, -0.1, 1.0).is_err()); // sigma<0
        assert!(StarlingParams::new(35.0, -2.0, 28.0, 5.0, 1.0, -1.0).is_err()); // Kf<0
        assert!(StarlingParams::new(f64::INFINITY, 0.0, 28.0, 5.0, 1.0, 1.0).is_err());
        assert!(StarlingParams::new(35.0, 0.0, f64::NAN, 5.0, 1.0, 1.0).is_err());
    }

    #[test]
    fn direction_tolerance_band_reads_as_equilibrium() {
        // A whisper of net pressure (+0.5) below a tol of 1.0 reads flat.
        let p = StarlingParams::new(23.5, 0.0, 28.0, 5.0, 1.0, 1.0).unwrap();
        assert!((p.net_driving_pressure() - 0.5).abs() < EPS);
        assert_eq!(p.direction_with_tol(1.0), FluxDirection::Equilibrium);
        // With a tight tolerance the same +0.5 reads as filtration.
        assert_eq!(p.direction_with_tol(0.0), FluxDirection::Filtration);
    }
}
