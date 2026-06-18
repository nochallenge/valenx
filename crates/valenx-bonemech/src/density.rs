//! Apparent-density strength scaling (a power law).
//!
//! Bone is porous, and both its stiffness and its strength rise steeply
//! with **apparent density** `rho` (the mass of bone tissue per unit of
//! *bulk* volume, including the pore space). The widely cited
//! Carter-Hayes family of relations models this as a power law: ultimate
//! strength scales roughly with the **square** of apparent density,
//!
//! ```text
//! sigma_ult(rho) = sigma_ref * (rho / rho_ref)^exponent
//! ```
//!
//! with `exponent ~ 2` as the textbook default. Denser bone is therefore
//! disproportionately stronger — halving the apparent density cuts the
//! predicted strength to about a quarter.
//!
//! ## Units
//!
//! | Quantity            | Symbol      | Unit    |
//! |---------------------|-------------|---------|
//! | Apparent density    | `rho`       | g/cm^3  |
//! | Ultimate stress     | `sigma_ult` | MPa     |
//! | Power-law exponent  | `exponent`  | —       |
//!
//! Density appears only as the **ratio** `rho / rho_ref`, so any
//! consistent mass-per-volume unit works as long as the queried and
//! reference densities share it; the unit cancels and the result inherits
//! the units of `sigma_ref` (MPa here).
//!
//! ## Model
//!
//! This is a single-power-law fit, not a micro-mechanical derivation.
//! Real exponents reported in the literature span roughly `2–3` and
//! differ between trabecular and cortical bone and between compression
//! and tension; the configurable [`PowerLaw::exponent`] lets a caller
//! match a specific published fit. Outside the calibrated density range
//! the extrapolation is unreliable.

use crate::error::{BoneError, Result};
use serde::{Deserialize, Serialize};

/// Default apparent-density exponent for the strength power law.
///
/// The textbook Carter-Hayes value: ultimate strength scales with the
/// square of apparent density.
pub const DEFAULT_DENSITY_EXPONENT: f64 = 2.0;

/// A calibrated apparent-density → strength power law.
///
/// Anchored at a single reference point `(rho_ref, sigma_ref)` with a
/// power-law `exponent`. Evaluate it with
/// [`strength`](PowerLaw::strength). Construct via [`PowerLaw::new`]
/// (explicit exponent) or [`PowerLaw::carter_hayes`] (the squared
/// default).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PowerLaw {
    /// Reference apparent density `rho_ref` (g/cm^3). Strictly positive.
    pub reference_density: f64,
    /// Ultimate stress `sigma_ref` (MPa) at the reference density.
    /// Strictly positive.
    pub reference_strength_mpa: f64,
    /// Power-law exponent. Non-negative and finite; the textbook value is
    /// [`DEFAULT_DENSITY_EXPONENT`] (2.0).
    pub exponent: f64,
}

impl PowerLaw {
    /// Build a power law from a reference point and an explicit exponent.
    ///
    /// # Errors
    ///
    /// Returns [`BoneError::Invalid`] if `reference_density <= 0`, if
    /// `reference_strength_mpa <= 0`, if `exponent < 0`, or if any
    /// argument is non-finite.
    pub fn new(reference_density: f64, reference_strength_mpa: f64, exponent: f64) -> Result<Self> {
        if !reference_density.is_finite() || reference_density <= 0.0 {
            return Err(BoneError::invalid(
                "reference_density",
                "reference apparent density must be positive and finite",
            ));
        }
        if !reference_strength_mpa.is_finite() || reference_strength_mpa <= 0.0 {
            return Err(BoneError::invalid(
                "reference_strength_mpa",
                "reference strength must be positive and finite",
            ));
        }
        if !exponent.is_finite() || exponent < 0.0 {
            return Err(BoneError::invalid(
                "exponent",
                "power-law exponent must be non-negative and finite",
            ));
        }
        Ok(PowerLaw {
            reference_density,
            reference_strength_mpa,
            exponent,
        })
    }

    /// Build a power law with the textbook squared exponent
    /// ([`DEFAULT_DENSITY_EXPONENT`] = 2.0) anchored at the given
    /// reference point.
    ///
    /// # Errors
    ///
    /// Returns [`BoneError::Invalid`] if `reference_density <= 0`,
    /// `reference_strength_mpa <= 0`, or either is non-finite.
    pub fn carter_hayes(reference_density: f64, reference_strength_mpa: f64) -> Result<Self> {
        PowerLaw::new(
            reference_density,
            reference_strength_mpa,
            DEFAULT_DENSITY_EXPONENT,
        )
    }

    /// Predicted ultimate stress at an apparent density `rho`, in MPa.
    ///
    /// Evaluates `sigma_ref * (rho / rho_ref)^exponent`.
    ///
    /// # Errors
    ///
    /// Returns [`BoneError::Invalid`] if `apparent_density < 0` or is
    /// non-finite. (A zero density yields zero strength for any positive
    /// exponent.)
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_bonemech::PowerLaw;
    /// // Reference 1.9 g/cm^3 -> 150 MPa, squared. Doubling the ratio to
    /// // 3.8 g/cm^3 multiplies strength by 2^2 = 4 -> 600 MPa.
    /// let law = PowerLaw::carter_hayes(1.9, 150.0).unwrap();
    /// let s = law.strength(3.8).unwrap();
    /// assert!((s - 600.0).abs() < 1e-9);
    /// ```
    pub fn strength(&self, apparent_density: f64) -> Result<f64> {
        if !apparent_density.is_finite() || apparent_density < 0.0 {
            return Err(BoneError::invalid(
                "apparent_density",
                "apparent density must be non-negative and finite",
            ));
        }
        let ratio = apparent_density / self.reference_density;
        Ok(self.reference_strength_mpa * ratio.powf(self.exponent))
    }

    /// Apparent density `rho` (in the reference-density units) that the law
    /// predicts for a target ultimate stress — the inverse of
    /// [`strength`](PowerLaw::strength).
    ///
    /// Inverting `sigma = sigma_ref * (rho / rho_ref)^exponent` gives
    ///
    /// ```text
    /// rho = rho_ref * (sigma / sigma_ref)^(1 / exponent)
    /// ```
    ///
    /// A zero target maps to zero density (for any positive exponent).
    ///
    /// # Errors
    ///
    /// Returns [`BoneError::Invalid`] if `target_strength_mpa` is negative
    /// or non-finite, or if the law's `exponent` is zero — a flat law has
    /// no density dependence, so the density cannot be recovered from the
    /// strength.
    pub fn density_for_strength(&self, target_strength_mpa: f64) -> Result<f64> {
        if !target_strength_mpa.is_finite() || target_strength_mpa < 0.0 {
            return Err(BoneError::invalid(
                "target_strength_mpa",
                "target ultimate stress must be non-negative and finite",
            ));
        }
        if self.exponent <= 0.0 {
            return Err(BoneError::invalid(
                "exponent",
                "strength is independent of density when the exponent is zero, \
                 so the density cannot be recovered",
            ));
        }
        let ratio = target_strength_mpa / self.reference_strength_mpa;
        Ok(self.reference_density * ratio.powf(1.0 / self.exponent))
    }
}

/// Predicted ultimate stress from apparent density via a power law, in
/// MPa.
///
/// A free-function shortcut around [`PowerLaw::strength`]:
///
/// ```text
/// sigma_ult = reference_strength_mpa
///           * (apparent_density / reference_density)^exponent
/// ```
///
/// Use `exponent = 2.0` ([`DEFAULT_DENSITY_EXPONENT`]) for the textbook
/// squared Carter-Hayes scaling. Denser bone is stronger: for any
/// positive exponent the result is strictly increasing in
/// `apparent_density`.
///
/// # Errors
///
/// Returns [`BoneError::Invalid`] if `reference_density <= 0`,
/// `reference_strength_mpa <= 0`, `exponent < 0`, `apparent_density < 0`,
/// or any argument is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_bonemech::strength_from_density;
/// // Halving the density quarters the strength under the squared law.
/// let full = strength_from_density(1.8, 1.8, 160.0, 2.0).unwrap();
/// let half = strength_from_density(0.9, 1.8, 160.0, 2.0).unwrap();
/// assert!((full - 160.0).abs() < 1e-9);
/// assert!((half - 40.0).abs() < 1e-9);
/// ```
pub fn strength_from_density(
    apparent_density: f64,
    reference_density: f64,
    reference_strength_mpa: f64,
    exponent: f64,
) -> Result<f64> {
    PowerLaw::new(reference_density, reference_strength_mpa, exponent)?.strength(apparent_density)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// At the reference density the law returns the reference strength
    /// exactly, for any exponent (the ratio is 1, and 1^k = 1).
    #[test]
    fn at_reference_density_returns_reference_strength() {
        for exp in [0.0, 1.0, 2.0, 2.5, 3.0] {
            let law = PowerLaw::new(1.9, 150.0, exp).unwrap();
            let s = law.strength(1.9).unwrap();
            assert!((s - 150.0).abs() < 1e-9, "exp {exp}: got {s}");
        }
    }

    /// The defining squared scaling: doubling apparent density multiplies
    /// strength by exactly 4.
    #[test]
    fn squared_law_doubling_quadruples_strength() {
        let law = PowerLaw::carter_hayes(1.0, 100.0).unwrap();
        let s1 = law.strength(1.0).unwrap();
        let s2 = law.strength(2.0).unwrap();
        assert!((s1 - 100.0).abs() < 1e-9, "got {s1}");
        assert!((s2 - 400.0).abs() < 1e-9, "got {s2}");
        assert!((s2 / s1 - 4.0).abs() < 1e-9);
    }

    /// Halving apparent density quarters strength under the squared law.
    #[test]
    fn squared_law_halving_quarters_strength() {
        let s = strength_from_density(0.9, 1.8, 160.0, 2.0).unwrap();
        assert!((s - 40.0).abs() < 1e-9, "got {s}");
    }

    /// Denser bone is stronger: the power law is strictly monotonically
    /// increasing across a sweep of densities (the headline property).
    #[test]
    fn denser_bone_is_strictly_stronger() {
        let law = PowerLaw::carter_hayes(1.5, 120.0).unwrap();
        let densities = [0.2, 0.5, 0.8, 1.1, 1.5, 1.9, 2.3];
        let mut prev = law.strength(densities[0]).unwrap();
        for &rho in &densities[1..] {
            let s = law.strength(rho).unwrap();
            assert!(
                s > prev,
                "strength not increasing at rho={rho}: {s} <= {prev}"
            );
            prev = s;
        }
    }

    /// The free function and the struct method agree exactly.
    #[test]
    fn free_function_matches_method() {
        let law = PowerLaw::new(1.7, 140.0, 2.3).unwrap();
        let via_method = law.strength(1.2).unwrap();
        let via_free = strength_from_density(1.2, 1.7, 140.0, 2.3).unwrap();
        assert!(
            (via_method - via_free).abs() < 1e-12,
            "{via_method} vs {via_free}"
        );
    }

    /// A cubic exponent scales as the cube of the density ratio.
    #[test]
    fn cubic_exponent_matches_analytic() {
        // ratio = 2 -> 2^3 = 8.
        let s = strength_from_density(2.0, 1.0, 50.0, 3.0).unwrap();
        assert!((s - 400.0).abs() < 1e-9, "got {s}");
    }

    /// Zero density yields zero strength under a positive exponent.
    #[test]
    fn zero_density_yields_zero_strength() {
        let s = strength_from_density(0.0, 1.8, 160.0, 2.0).unwrap();
        assert!(s.abs() < 1e-12, "got {s}");
    }

    #[test]
    fn invalid_inputs_rejected() {
        assert!(PowerLaw::new(0.0, 150.0, 2.0).is_err());
        assert!(PowerLaw::new(1.9, 0.0, 2.0).is_err());
        assert!(PowerLaw::new(1.9, 150.0, -1.0).is_err());
        assert!(PowerLaw::new(f64::NAN, 150.0, 2.0).is_err());
        let law = PowerLaw::carter_hayes(1.9, 150.0).unwrap();
        assert!(law.strength(-0.1).is_err());
        assert!(law.strength(f64::INFINITY).is_err());
        assert!(strength_from_density(1.0, -1.0, 150.0, 2.0).is_err());
    }

    /// The inverse recovers the density that produces a target strength,
    /// in both directions, for the squared law.
    #[test]
    fn density_for_strength_inverts_strength() {
        let law = PowerLaw::carter_hayes(1.9, 150.0).unwrap();
        // Forward then inverse returns the original density.
        for &rho in &[0.5, 1.0, 1.9, 2.4] {
            let sigma = law.strength(rho).unwrap();
            let back = law.density_for_strength(sigma).unwrap();
            assert!((back - rho).abs() < 1e-9, "rho {rho} -> {back}");
        }
        // Inverse then forward returns the original strength.
        for &sigma in &[40.0, 150.0, 600.0] {
            let rho = law.density_for_strength(sigma).unwrap();
            let fwd = law.strength(rho).unwrap();
            assert!((fwd - sigma).abs() < 1e-7, "sigma {sigma} -> {fwd}");
        }
    }

    /// Hand value and the reference-point fixed point.
    #[test]
    fn density_for_strength_matches_hand_value() {
        // Squared law anchored at (1.0, 100.0): sigma 400 -> rho
        // 1 * (400/100)^(1/2) = 2.
        let law = PowerLaw::carter_hayes(1.0, 100.0).unwrap();
        let rho = law.density_for_strength(400.0).unwrap();
        assert!((rho - 2.0).abs() < 1e-9, "got {rho}");
        // At the reference strength the inverse returns the reference density.
        let at_ref = law.density_for_strength(100.0).unwrap();
        assert!((at_ref - 1.0).abs() < 1e-9, "got {at_ref}");
    }

    /// Zero target maps to zero density; bad targets and a flat (zero
    /// exponent) law are rejected.
    #[test]
    fn density_for_strength_handles_zero_and_rejects_bad_input() {
        let law = PowerLaw::carter_hayes(1.8, 160.0).unwrap();
        // Zero target strength -> zero density (for a positive exponent).
        assert!(law.density_for_strength(0.0).unwrap().abs() < 1e-12);
        // Negative or non-finite target is rejected.
        assert!(law.density_for_strength(-1.0).is_err());
        assert!(law.density_for_strength(f64::NAN).is_err());
        // A zero-exponent law has no density dependence, so the inverse is
        // undefined.
        let flat = PowerLaw::new(1.8, 160.0, 0.0).unwrap();
        assert!(flat.density_for_strength(160.0).is_err());
    }

    /// The inverse works for a non-squared (cubic) exponent.
    #[test]
    fn density_for_strength_works_for_cubic_exponent() {
        // Cubic law: ratio 2 -> 8x strength, so inverting 8x recovers
        // ratio 2.
        let law = PowerLaw::new(1.0, 50.0, 3.0).unwrap();
        let rho = law.density_for_strength(400.0).unwrap();
        assert!((rho - 2.0).abs() < 1e-9, "got {rho}");
    }
}
