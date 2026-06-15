//! Buoyancy — Archimedes' principle and floating equilibrium.
//!
//! A body immersed in a fluid feels an upward buoyant force equal to the
//! weight of the fluid it displaces (Archimedes' principle):
//!
//! ```text
//! F_b = rho_fluid * g * V_displaced
//! ```
//!
//! For a *fully submerged* body of volume `V`, `V_displaced = V`. A
//! *freely floating* body sinks until the displaced weight balances its
//! own weight, so it displaces the fraction
//!
//! ```text
//! f = rho_body / rho_fluid          (the submerged / floating fraction)
//! ```
//!
//! of its own volume. The body floats when `f < 1` (i.e. it is less
//! dense than the fluid), is neutrally buoyant when `f == 1`, and sinks
//! when `f > 1`. The classic "≈ 90 % of an iceberg is underwater"
//! follows directly: `917 / 1025 ≈ 0.895`.

use crate::error::{require_non_negative, require_positive, FluidStaticsError, Result};
use crate::fluid::{Fluid, STANDARD_GRAVITY};

/// The buoyant (upward) force on a body that displaces `volume_disp_m3`
/// of `fluid` under gravity `gravity`, in newtons:
/// `F_b = rho_fluid * g * V_displaced`.
///
/// This is the force regardless of whether the body is fully submerged
/// or floating — only the *displaced* volume matters.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gravity` is not strictly positive or `volume_disp_m3` is negative
/// / non-finite.
pub fn buoyant_force(fluid: &Fluid, gravity: f64, volume_disp_m3: f64) -> Result<f64> {
    let gravity = require_positive("gravity", gravity)?;
    let volume = require_non_negative("volume_disp_m3", volume_disp_m3)?;
    Ok(fluid.density() * gravity * volume)
}

/// The weight of a body of `mass_kg` under gravity `gravity`, in newtons:
/// `W = m * g`.
///
/// Provided so callers can compare a body's weight against its buoyant
/// force without reaching for the bare formula.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gravity` is not strictly positive or `mass_kg` is negative /
/// non-finite.
pub fn weight(mass_kg: f64, gravity: f64) -> Result<f64> {
    let gravity = require_positive("gravity", gravity)?;
    let mass = require_non_negative("mass_kg", mass_kg)?;
    Ok(mass * gravity)
}

/// The net upward force on a fully submerged body, in newtons:
/// `F_net = F_buoyant - W = (rho_fluid - rho_body) * g * V`.
///
/// A positive value means the body rises (it is less dense than the
/// fluid), zero means it is neutrally buoyant, and a negative value
/// means it sinks.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gravity` is not strictly positive or any density / volume is out
/// of its physical domain.
pub fn net_force_submerged(
    body: &Fluid,
    fluid: &Fluid,
    gravity: f64,
    volume_m3: f64,
) -> Result<f64> {
    let gravity = require_positive("gravity", gravity)?;
    let volume = require_non_negative("volume_m3", volume_m3)?;
    Ok((fluid.density() - body.density()) * gravity * volume)
}

/// Outcome of placing a body of density `rho_body` into a fluid of
/// density `rho_fluid`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum FloatState {
    /// The body is less dense than the fluid and floats partially
    /// submerged.
    Floats,
    /// The body has exactly the fluid's density and is in neutral
    /// equilibrium at any depth.
    Neutral,
    /// The body is denser than the fluid and sinks.
    Sinks,
}

/// Result of a free-floating equilibrium analysis for a homogeneous body
/// of known density and volume placed in a fluid.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FloatResult {
    /// Whether the body floats, is neutral, or sinks.
    pub state: FloatState,
    /// The fraction of the body's volume that lies below the free
    /// surface at equilibrium, `rho_body / rho_fluid`, clamped to the
    /// physically meaningful range `[0, 1]` for a floating body. For a
    /// sinking body the *unclamped* ratio exceeds 1; the clamped value
    /// is reported as 1 (the body is fully submerged).
    pub submerged_fraction: f64,
    /// The displaced volume at equilibrium, in cubic metres. For a
    /// floating body this is `submerged_fraction * V`; for a sinking or
    /// neutral body it is the full body volume `V`.
    pub displaced_volume_m3: f64,
}

impl FloatResult {
    /// The fraction of the body's volume that sits *above* the free
    /// surface (the visible freeboard), `1 - submerged_fraction`. Zero
    /// for a sinking or neutrally buoyant body.
    pub fn exposed_fraction(&self) -> f64 {
        1.0 - self.submerged_fraction
    }

    /// Whether the body floats (any visible freeboard at all).
    pub fn is_floating(&self) -> bool {
        self.state == FloatState::Floats
    }
}

/// The submerged fraction of a freely floating homogeneous body,
/// `rho_body / rho_fluid` — the share of its volume that lies below the
/// waterline at equilibrium.
///
/// Equilibrium requires the displaced weight to equal the body weight:
/// `rho_fluid * g * (f * V) = rho_body * g * V`, hence
/// `f = rho_body / rho_fluid`. The raw ratio is returned *unclamped*, so
/// a value greater than 1 signals a body that cannot float (it would
/// have to displace more than its own volume).
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if either density is out of its physical domain (densities come from
/// [`Fluid`], so this only guards against a future non-positive source).
pub fn floating_fraction(body: &Fluid, fluid: &Fluid) -> Result<f64> {
    // `Fluid` guarantees positive density, but guard defensively so the
    // function is correct for any density source.
    let body_rho = require_positive("body_density", body.density())?;
    let fluid_rho = require_positive("fluid_density", fluid.density())?;
    Ok(body_rho / fluid_rho)
}

/// Analyse the free-floating equilibrium of a homogeneous body of the
/// given density and `volume_m3` placed in `fluid`.
///
/// Returns whether it floats / is neutral / sinks, the submerged
/// fraction (`rho_body / rho_fluid`, clamped to `[0, 1]`), and the
/// displaced volume at equilibrium.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `volume_m3` is negative / non-finite or a density is out of domain.
pub fn analyse_float(body: &Fluid, fluid: &Fluid, volume_m3: f64) -> Result<FloatResult> {
    let volume = require_non_negative("volume_m3", volume_m3)?;
    let raw = floating_fraction(body, fluid)?;

    // Compare densities with a relative tolerance so a body at exactly
    // the fluid density is reported as neutral rather than as a
    // razor-thin float or sink.
    let rel = (body.density() - fluid.density()).abs() / fluid.density();
    let state = if rel <= 1e-12 {
        FloatState::Neutral
    } else if raw < 1.0 {
        FloatState::Floats
    } else {
        FloatState::Sinks
    };

    let submerged_fraction = match state {
        FloatState::Floats => raw,
        FloatState::Neutral | FloatState::Sinks => 1.0,
    };
    let displaced_volume_m3 = submerged_fraction * volume;

    Ok(FloatResult {
        state,
        submerged_fraction,
        displaced_volume_m3,
    })
}

/// The waterline depth (draft) of a homogeneous body of uniform
/// horizontal cross-section `waterplane_area_m2` and total `volume_m3`
/// floating in `fluid`, in metres.
///
/// For a prism of constant cross-section the submerged depth is
/// `draft = submerged_fraction * V / A = f * height`. Useful for
/// reading how deep a uniform block or pontoon sits.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `waterplane_area_m2` is not strictly positive or `volume_m3` is
/// negative / non-finite, and
/// [`Geometry`](crate::FluidStaticsError::Geometry)
/// if the body does not float (a draft is only defined for a floater).
pub fn floating_draft(
    body: &Fluid,
    fluid: &Fluid,
    volume_m3: f64,
    waterplane_area_m2: f64,
) -> Result<f64> {
    let area = require_positive("waterplane_area_m2", waterplane_area_m2)?;
    let result = analyse_float(body, fluid, volume_m3)?;
    if result.state == FloatState::Sinks {
        return Err(FluidStaticsError::geometry(
            "floating draft",
            "body sinks (denser than fluid); draft is undefined",
        ));
    }
    Ok(result.displaced_volume_m3 / area)
}

/// Buoyant force on a body that displaces `volume_disp_m3` of `fluid`
/// under [`STANDARD_GRAVITY`], in newtons — a convenience wrapper around
/// [`buoyant_force`].
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `volume_disp_m3` is negative / non-finite.
pub fn buoyant_force_standard(fluid: &Fluid, volume_disp_m3: f64) -> Result<f64> {
    buoyant_force(fluid, STANDARD_GRAVITY, volume_disp_m3)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn buoyant_force_equals_weight_of_displaced_fluid() {
        // 2 m^3 displaced in water -> weight of 2000 kg of water.
        let f = Fluid::water();
        let fb = buoyant_force(&f, STANDARD_GRAVITY, 2.0).unwrap();
        let displaced_mass = f.density() * 2.0; // 2000 kg
        let displaced_weight = weight(displaced_mass, STANDARD_GRAVITY).unwrap();
        assert!((fb - displaced_weight).abs() < 1e-6, "fb={fb}");
    }

    #[test]
    fn floating_body_displaces_its_own_weight() {
        // A floating body's buoyant force must exactly equal its weight.
        let body = Fluid::new(600.0).unwrap(); // e.g. light wood
        let fluid = Fluid::water();
        let volume = 0.5_f64; // m^3
        let body_mass = body.density() * volume;
        let body_weight = weight(body_mass, STANDARD_GRAVITY).unwrap();

        let res = analyse_float(&body, &fluid, volume).unwrap();
        let fb = buoyant_force(&fluid, STANDARD_GRAVITY, res.displaced_volume_m3).unwrap();
        assert!((fb - body_weight).abs() < 1e-6, "fb={fb} W={body_weight}");
    }

    #[test]
    fn floating_fraction_is_density_ratio() {
        let body = Fluid::new(800.0).unwrap();
        let fluid = Fluid::water();
        let f = floating_fraction(&body, &fluid).unwrap();
        assert!((f - 0.8).abs() < EPS, "got {f}");
    }

    #[test]
    fn floats_when_less_dense_than_fluid() {
        let body = Fluid::new(700.0).unwrap();
        let res = analyse_float(&body, &Fluid::water(), 1.0).unwrap();
        assert_eq!(res.state, FloatState::Floats);
        assert!(res.is_floating());
        assert!((res.submerged_fraction - 0.7).abs() < EPS);
        assert!((res.exposed_fraction() - 0.3).abs() < EPS);
    }

    #[test]
    fn sinks_when_denser_than_fluid() {
        let body = Fluid::new(2700.0).unwrap(); // aluminium in water
        let res = analyse_float(&body, &Fluid::water(), 1.0).unwrap();
        assert_eq!(res.state, FloatState::Sinks);
        assert!(!res.is_floating());
        // Sinking body is fully submerged: clamped fraction is 1.
        assert!((res.submerged_fraction - 1.0).abs() < EPS);
        // Unclamped ratio is well above 1.
        let raw = floating_fraction(&body, &Fluid::water()).unwrap();
        assert!(raw > 2.0, "got {raw}");
    }

    #[test]
    fn neutral_when_equal_density() {
        let body = Fluid::water();
        let res = analyse_float(&body, &Fluid::water(), 1.0).unwrap();
        assert_eq!(res.state, FloatState::Neutral);
        assert!((res.submerged_fraction - 1.0).abs() < EPS);
    }

    #[test]
    fn net_force_sign_tracks_density_difference() {
        // Less dense than fluid -> positive (rises).
        let up = net_force_submerged(
            &Fluid::new(500.0).unwrap(),
            &Fluid::water(),
            STANDARD_GRAVITY,
            1.0,
        )
        .unwrap();
        assert!(up > 0.0, "got {up}");

        // Denser than fluid -> negative (sinks).
        let down = net_force_submerged(
            &Fluid::new(1500.0).unwrap(),
            &Fluid::water(),
            STANDARD_GRAVITY,
            1.0,
        )
        .unwrap();
        assert!(down < 0.0, "got {down}");

        // Equal density -> zero.
        let zero =
            net_force_submerged(&Fluid::water(), &Fluid::water(), STANDARD_GRAVITY, 1.0).unwrap();
        assert!(zero.abs() < 1e-6, "got {zero}");
    }

    #[test]
    fn iceberg_is_about_ninety_percent_submerged() {
        // Sea ice ~917 kg/m^3 in seawater ~1025 kg/m^3.
        let ice = Fluid::new(917.0).unwrap();
        let res = analyse_float(&ice, &Fluid::seawater(), 1.0).unwrap();
        assert_eq!(res.state, FloatState::Floats);
        assert!(
            (res.submerged_fraction - 0.8946).abs() < 1e-3,
            "got {}",
            res.submerged_fraction
        );
    }

    #[test]
    fn floating_draft_is_fraction_times_height() {
        // A uniform prism: V = A * H, so draft = f * H.
        let body = Fluid::new(600.0).unwrap();
        let fluid = Fluid::water();
        let area = 2.0_f64; // m^2
        let height = 0.5_f64; // m
        let volume = area * height;
        let draft = floating_draft(&body, &fluid, volume, area).unwrap();
        // f = 0.6 -> draft = 0.3 m.
        assert!((draft - 0.3).abs() < EPS, "got {draft}");
    }

    #[test]
    fn floating_draft_errors_for_sinking_body() {
        let body = Fluid::new(2000.0).unwrap();
        let err = floating_draft(&body, &Fluid::water(), 1.0, 1.0);
        assert!(err.is_err());
    }

    #[test]
    fn rejects_bad_arguments() {
        assert!(buoyant_force(&Fluid::water(), -1.0, 1.0).is_err());
        assert!(buoyant_force(&Fluid::water(), STANDARD_GRAVITY, -1.0).is_err());
        assert!(weight(-1.0, STANDARD_GRAVITY).is_err());
        assert!(floating_draft(&Fluid::new(500.0).unwrap(), &Fluid::water(), 1.0, 0.0).is_err());
    }

    #[test]
    fn standard_gravity_wrapper_matches_explicit() {
        let a = buoyant_force_standard(&Fluid::water(), 1.5).unwrap();
        let b = buoyant_force(&Fluid::water(), STANDARD_GRAVITY, 1.5).unwrap();
        assert!((a - b).abs() < 1e-9);
    }
}
