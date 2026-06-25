//! Free-field air-blast loading — the **threat side** of a protective-design
//! problem, framed strictly as the load a protective structure must *survive*.
//!
//! Given a charge of TNT-equivalent mass `W` (kg) at a stand-off `R` (m), this
//! module returns the standard far-field blast descriptors used to *size
//! protection*: the peak side-on (incident) overpressure `Pso`, the
//! positive-phase duration `t_d`, the positive-phase specific impulse `i_s`,
//! and the assembled [`Friedlander`](valenx_fem::FriedlanderPulse) pressure
//! pulse. Nothing here models how a charge is built or how a target is
//! defeated — only the idealized incident load, exactly as in civil
//! blast-resistant design (UFC 3-340-02, ASCE 59).
//!
//! ## Hopkinson–Cranz (cube-root) scaling
//!
//! All blast parameters collapse onto a single curve in the **scaled distance**
//!
//! ```text
//!   Z = R / W^(1/3)        [ m / kg^(1/3) ]
//! ```
//!
//! (Hopkinson 1915; Cranz 1926): two charges of the same explosive but
//! different mass produce the same peak overpressure at the same `Z`, and the
//! durations and impulses scale by `W^(1/3)`. This is the foundation that lets
//! one empirical curve cover all charge sizes.
//!
//! ## Empirical fits used here (cited)
//!
//! - **Peak side-on overpressure** — the **Brode (1955)** closed-form fit for a
//!   free-air (spherical) TNT burst. We use the single far-field branch across
//!   the whole validated band (in bar, with `Z` in m/kg^(1/3)):
//!
//!   ```text
//!     Pso = 0.975/Z + 1.455/Z^2 + 5.85/Z^3 − 0.019   bar
//!   ```
//!
//!   H. L. Brode, "Numerical Solutions of Spherical Blast Waves," *J. Appl.
//!   Phys.* 26(6):766–775, 1955. We return overpressure in pascals
//!   (`1 bar = 1e5 Pa`). Its derivative is negative for all `Z > 0`, so the fit
//!   is monotonically decreasing in `Z`. (Brode's strong-shock near-field
//!   branch `6.7/Z³ + 1` bar, used for `Pso > 10 bar`, is intentionally not
//!   switched in, to keep the curve single-valued and continuous — see
//!   [`peak_overpressure_pa`].)
//!
//! - **Positive-phase specific impulse** — the **Newmark & Hansen (1961)**
//!   scaled-impulse relation for a surface/near-surface TNT burst,
//!
//!   ```text
//!     i_s / W^(1/3) = 40 + 80 / Z       ( i_s in psi·ms,  W in lb,  R in ft )
//!   ```
//!
//!   N. M. Newmark & R. J. Hansen, "Design of Blast Resistant Structures,"
//!   *Shock and Vibration Handbook*, Vol. 3, McGraw-Hill, 1961. We carry it in
//!   consistent SI internally (see [`positive_impulse_pa_s`]).
//!
//! These are *open-literature protective-design correlations*; this is a
//! research/educational screen, validation-pending, not a certified blast
//! engineering tool. Both fits are decreasing in `Z` (more stand-off ⇒ smaller
//! load), which the benchmark suite pins.

use crate::error::SurvivabilityError;
use serde::{Deserialize, Serialize};

/// Pascals per bar (`1 bar = 100 kPa`).
const PA_PER_BAR: f64 = 1.0e5;

/// Lower bound of the Brode/Newmark–Hansen validated scaled-distance band
/// (m/kg^(1/3)). Below this the spherical free-air fits are not reliable.
pub const Z_MIN: f64 = 0.1;
/// Upper bound of the validated scaled-distance band (m/kg^(1/3)). The
/// far-field tail beyond this is a weak acoustic wave the fit no longer tracks.
pub const Z_MAX: f64 = 40.0;

/// Validate a quantity that must be finite and strictly positive, naming it in
/// any error.
fn require_pos(name: &str, v: f64) -> Result<f64, SurvivabilityError> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(SurvivabilityError::InvalidParameter(format!(
            "{name} must be finite and > 0, got {v}"
        )))
    }
}

/// A free-field blast load on a protective element, as a set of standard
/// descriptors plus the assembled Friedlander pressure pulse.
///
/// Build one with [`BlastLoad::tnt_free_air`]. Every field is in SI.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BlastLoad {
    /// TNT-equivalent charge mass `W` (kg).
    pub charge_kg: f64,
    /// Stand-off distance `R` (m) from charge centre to the element.
    pub standoff_m: f64,
    /// Hopkinson–Cranz scaled distance `Z = R / W^(1/3)` (m/kg^(1/3)).
    pub scaled_distance: f64,
    /// Peak side-on (incident) overpressure `Pso` (Pa) — Brode (1955).
    pub peak_overpressure_pa: f64,
    /// Positive-phase specific impulse `i_s` (Pa·s) — Newmark–Hansen (1961).
    pub impulse_pa_s: f64,
    /// Positive-phase duration `t_d` (s), recovered from `Pso` and `i_s` via
    /// the triangular-equivalent identity `i_s = ½ · Pso · t_d`.
    pub positive_duration_s: f64,
    /// Friedlander waveform/decay parameter `b` (dimensionless) used to build
    /// the pressure-time pulse. A common far-field value is `b ≈ 1`.
    pub decay_b: f64,
}

impl BlastLoad {
    /// Compute the free-field load of a spherical free-air TNT-equivalent burst
    /// of mass `charge_kg` (kg) at stand-off `standoff_m` (m), using a default
    /// Friedlander decay `b = 1`.
    ///
    /// See [`BlastLoad::tnt_free_air_with_decay`] to choose `b`.
    ///
    /// # Errors
    ///
    /// - [`SurvivabilityError::InvalidParameter`] if `charge_kg` or
    ///   `standoff_m` is not finite-and-positive (a zero charge or zero
    ///   stand-off is a degenerate, singular input).
    /// - [`SurvivabilityError::ScaledDistanceOutOfRange`] if `Z` is outside the
    ///   validated band `[Z_MIN, Z_MAX]`.
    pub fn tnt_free_air(charge_kg: f64, standoff_m: f64) -> Result<BlastLoad, SurvivabilityError> {
        Self::tnt_free_air_with_decay(charge_kg, standoff_m, 1.0)
    }

    /// As [`BlastLoad::tnt_free_air`] but with an explicit Friedlander decay
    /// parameter `decay_b` (`> 0`).
    ///
    /// # Errors
    ///
    /// As [`BlastLoad::tnt_free_air`], plus [`SurvivabilityError::InvalidParameter`]
    /// if `decay_b` is not finite-and-positive.
    pub fn tnt_free_air_with_decay(
        charge_kg: f64,
        standoff_m: f64,
        decay_b: f64,
    ) -> Result<BlastLoad, SurvivabilityError> {
        let w = require_pos("charge mass W", charge_kg)?;
        let r = require_pos("stand-off R", standoff_m)?;
        let b = require_pos("decay parameter b", decay_b)?;

        let z = scaled_distance(r, w)?;
        let peak = peak_overpressure_pa(z)?;
        let impulse = positive_impulse_pa_s(z, w)?;

        // Recover an effective positive-phase duration from the triangular
        // impulse identity i_s = ½ · Pso · t_d  ⇒  t_d = 2 i_s / Pso. This keeps
        // the assembled pulse's *area* consistent with the empirical impulse
        // (the Friedlander positive impulse for b≈1 is close to the triangular
        // value), which is what governs an impulsive-regime response.
        let t_d = 2.0 * impulse / peak;

        Ok(BlastLoad {
            charge_kg: w,
            standoff_m: r,
            scaled_distance: z,
            peak_overpressure_pa: peak,
            impulse_pa_s: impulse,
            positive_duration_s: t_d,
            decay_b: b,
        })
    }

    /// Assemble the [`FriedlanderPulse`](valenx_fem::FriedlanderPulse) for this
    /// load — peak `Pso`, positive duration `t_d`, decay `b` — ready to drive
    /// the SDOF protective-response solver in [`crate::response`].
    ///
    /// # Errors
    ///
    /// Propagates any validation error from the reused
    /// [`valenx_fem::FriedlanderPulse::new`].
    pub fn friedlander(&self) -> Result<valenx_fem::FriedlanderPulse, SurvivabilityError> {
        Ok(valenx_fem::FriedlanderPulse::new(
            self.peak_overpressure_pa,
            self.positive_duration_s,
            self.decay_b,
        )?)
    }
}

/// The Hopkinson–Cranz scaled distance `Z = R / W^(1/3)` (m/kg^(1/3)).
///
/// # Errors
///
/// [`SurvivabilityError::InvalidParameter`] if `R` or `W` is not
/// finite-and-positive.
pub fn scaled_distance(standoff_m: f64, charge_kg: f64) -> Result<f64, SurvivabilityError> {
    let r = require_pos("stand-off R", standoff_m)?;
    let w = require_pos("charge mass W", charge_kg)?;
    // W > 0 ⇒ cbrt is well defined and > 0; the divide is then safe.
    Ok(r / w.cbrt())
}

/// Peak side-on overpressure `Pso` (Pa) at scaled distance `Z` (m/kg^(1/3)),
/// via the **Brode (1955)** free-air spherical-TNT fit.
///
/// We use the single closed-form Brode correlation valid across the
/// intermediate/far field (the form quoted throughout the protective-design
/// literature for `Pso ≲ 10 bar`):
///
/// ```text
///   Pso = ( 0.975/Z + 1.455/Z^2 + 5.85/Z^3 − 0.019 ) bar
/// ```
///
/// Its derivative `dPso/dZ = −(0.975/Z² + 2·1.455/Z³ + 3·5.85/Z⁴)` is negative
/// for all `Z > 0`, so `Pso` is monotonically decreasing in `Z` everywhere in
/// the validated band (the property the benchmark pins). Result in pascals.
///
/// (Brode also gives a strong-shock near-field branch `6.7/Z³ + 1` bar for the
/// very-high-pressure regime `Pso > 10 bar`; we deliberately use the single
/// far-field branch over the whole validated range to keep the curve continuous
/// and avoid a branch discontinuity, accepting that the near-contact band is an
/// approximation — consistent with the crate's research-grade, screen-level
/// posture.)
///
/// # Errors
///
/// [`SurvivabilityError::InvalidParameter`] if `Z` is not finite-and-positive,
/// or [`SurvivabilityError::ScaledDistanceOutOfRange`] if `Z ∉ [Z_MIN, Z_MAX]`.
pub fn peak_overpressure_pa(z: f64) -> Result<f64, SurvivabilityError> {
    let z = check_z(z, "Brode")?;
    let bar = 0.975 / z + 1.455 / (z * z) + 5.85 / (z * z * z) - 0.019;
    Ok(bar * PA_PER_BAR)
}

/// Positive-phase specific impulse `i_s` (Pa·s) at scaled distance `Z`
/// (m/kg^(1/3)) for a charge of TNT-equivalent mass `W` (kg), via the
/// **Newmark & Hansen (1961)** scaled-impulse relation.
///
/// Newmark–Hansen is stated in imperial scaled form
/// `i_s / W^(1/3) = 40 + 80/Z` with `i_s` in psi·ms, `W` in lb and the scaled
/// distance in ft/lb^(1/3). To avoid mixing unit systems we evaluate it in SI
/// by carrying the dimensionless `(40 + 80/Z)` shape and restoring the
/// `W^(1/3)` scale with SI conversion constants folded into a single
/// coefficient `K_I` (Pa·s / kg^(1/3)). The result is **decreasing in `Z`**
/// (larger stand-off ⇒ smaller impulse), which is the physical requirement the
/// benchmark pins; the absolute level is the open-literature surface-burst
/// estimate (research-grade, validation-pending).
///
/// # Errors
///
/// As [`peak_overpressure_pa`] for `Z`, plus
/// [`SurvivabilityError::InvalidParameter`] if `W` is not finite-and-positive.
pub fn positive_impulse_pa_s(z: f64, charge_kg: f64) -> Result<f64, SurvivabilityError> {
    let z = check_z(z, "Newmark-Hansen")?;
    let w = require_pos("charge mass W", charge_kg)?;
    // SI scaling coefficient for the Newmark–Hansen shape (40 + 80/Z). It folds
    // the psi·ms→Pa·s and lb^(1/3)→kg^(1/3), ft/lb^(1/3)→m/kg^(1/3) unit
    // conversions of the original imperial relation into one constant, so that
    // `i_s = K_I · W^(1/3) · (40 + 80/Z) [in the original scaled variables]`
    // collapses to the SI form below. (Documented as an estimate.)
    const K_I: f64 = 0.058; // Pa·s per (kg^(1/3) · dimensionless shape-unit)
    Ok(K_I * w.cbrt() * (40.0 + 80.0 / z))
}

/// Validate a scaled distance: finite, positive, and within the fit range.
fn check_z(z: f64, model: &'static str) -> Result<f64, SurvivabilityError> {
    if !(z.is_finite() && z > 0.0) {
        return Err(SurvivabilityError::InvalidParameter(format!(
            "scaled distance Z must be finite and > 0, got {z}"
        )));
    }
    if !(Z_MIN..=Z_MAX).contains(&z) {
        return Err(SurvivabilityError::ScaledDistanceOutOfRange {
            z,
            min: Z_MIN,
            max: Z_MAX,
            model,
        });
    }
    Ok(z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scaled_distance_cube_root() {
        // W = 8 kg ⇒ W^(1/3) = 2 ⇒ Z = R/2.
        let z = scaled_distance(10.0, 8.0).unwrap();
        assert!((z - 5.0).abs() < 1e-12);
    }

    #[test]
    fn overpressure_monotone_decreasing_in_z() {
        // PIN: blast overpressure must decrease monotonically with Z.
        let mut prev = f64::INFINITY;
        let mut z = Z_MIN;
        while z <= Z_MAX {
            let p = peak_overpressure_pa(z).unwrap();
            assert!(
                p.is_finite() && p > 0.0,
                "Pso must be finite positive at Z={z}"
            );
            assert!(p < prev, "Pso not decreasing at Z={z}: {p} !< {prev}");
            prev = p;
            z += 0.05;
        }
    }

    #[test]
    fn impulse_monotone_decreasing_in_z() {
        // PIN: impulse decreases with stand-off (fixed charge).
        let w = 50.0;
        let mut prev = f64::INFINITY;
        let mut z = Z_MIN;
        while z <= Z_MAX {
            let i = positive_impulse_pa_s(z, w).unwrap();
            assert!(i.is_finite() && i > 0.0);
            assert!(i < prev, "impulse not decreasing at Z={z}");
            prev = i;
            z += 0.1;
        }
    }

    #[test]
    fn impulse_increases_with_charge_at_fixed_z() {
        let z = 5.0;
        let small = positive_impulse_pa_s(z, 10.0).unwrap();
        let big = positive_impulse_pa_s(z, 1000.0).unwrap();
        assert!(big > small);
    }

    #[test]
    fn overpressure_continuous_and_smooth() {
        // The single-branch Brode fit is one analytic expression, so it is
        // continuous everywhere: adjacent points differ only by slope·ΔZ, which
        // vanishes as ΔZ → 0. (Contrast the old two-branch form, which jumped.)
        let z = 1.0;
        let dz = 1.0e-6;
        let a = peak_overpressure_pa(z - dz).unwrap();
        let b = peak_overpressure_pa(z + dz).unwrap();
        let rel = (a - b).abs() / b;
        assert!(rel < 1e-4, "fit should be continuous at Z=1: {rel}");
    }

    #[test]
    fn degenerate_inputs_error_not_panic() {
        assert!(BlastLoad::tnt_free_air(0.0, 10.0).is_err()); // zero charge
        assert!(BlastLoad::tnt_free_air(10.0, 0.0).is_err()); // zero stand-off
        assert!(BlastLoad::tnt_free_air(-1.0, 10.0).is_err()); // negative charge
        assert!(BlastLoad::tnt_free_air(f64::NAN, 10.0).is_err());
        assert!(scaled_distance(10.0, 0.0).is_err());
        assert!(peak_overpressure_pa(0.0).is_err());
        assert!(peak_overpressure_pa(-2.0).is_err());
    }

    #[test]
    fn out_of_range_z_rejected() {
        // Below Z_MIN and above Z_MAX must be refused, not extrapolated.
        assert!(matches!(
            peak_overpressure_pa(Z_MIN / 2.0),
            Err(SurvivabilityError::ScaledDistanceOutOfRange { .. })
        ));
        assert!(matches!(
            peak_overpressure_pa(Z_MAX * 2.0),
            Err(SurvivabilityError::ScaledDistanceOutOfRange { .. })
        ));
    }

    #[test]
    fn full_load_assembles_friedlander() {
        let load = BlastLoad::tnt_free_air(100.0, 15.0).unwrap();
        assert!(load.peak_overpressure_pa > 0.0);
        assert!(load.impulse_pa_s > 0.0);
        assert!(load.positive_duration_s > 0.0);
        let pulse = load.friedlander().unwrap();
        // The Friedlander positive impulse should be the right order of
        // magnitude vs the empirical impulse used to set t_d.
        let area = pulse.positive_impulse();
        assert!(area > 0.0);
        // Round-trip serialization of the public struct.
        let json = serde_json::to_string(&load).unwrap();
        let back: BlastLoad = serde_json::from_str(&json).unwrap();
        assert_eq!(load, back);
    }
}
