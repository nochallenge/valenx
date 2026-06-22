//! Black-hole observables: horizons, ergosphere, photon sphere, ISCO, shadow
//! radius and gravitational redshift.
//!
//! Geometrized units (`G = c = 1`) throughout. Each quantity uses a closed form
//! where one exists; for parameter combinations without an implemented closed
//! form (e.g. an ISCO with *both* spin and charge) the functions return
//! [`RelativityError::Unsupported`] rather than guessing — the numerical
//! geodesic path (a later layer) covers those.

use std::f64::consts::FRAC_PI_2;

use crate::metric::Spacetime;
use crate::spacetimes::KerrNewman;
use crate::{RelativityError, Result};

/// Sense of an equatorial orbit relative to the hole's spin.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum OrbitSense {
    /// Co-rotating with the black hole's spin.
    Prograde,
    /// Counter-rotating against the spin.
    Retrograde,
}

impl OrbitSense {
    /// `+1` for prograde, `−1` for retrograde — the sign multiplying `a` in the
    /// equatorial orbit formulas.
    fn sign(self) -> f64 {
        match self {
            OrbitSense::Prograde => 1.0,
            OrbitSense::Retrograde => -1.0,
        }
    }
}

/// Event- and inner-horizon radii of a black hole.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Horizons {
    /// Outer (event) horizon `r+ = M + √(M² − a² − Q²)`.
    pub outer: f64,
    /// Inner (Cauchy) horizon `r− = M − √(M² − a² − Q²)`.
    pub inner: f64,
}

/// Validate that `bh` is a genuine (sub-extremal, positive-mass) black hole.
fn check_blackhole(bh: &KerrNewman) -> Result<()> {
    if !bh.mass.is_finite() || bh.mass <= 0.0 {
        return Err(RelativityError::InvalidParameter(format!(
            "mass must be a positive finite number, got {}",
            bh.mass
        )));
    }
    if !bh.is_subextremal() {
        return Err(RelativityError::InvalidParameter(
            "super-extremal a²+Q² > M²: a naked singularity, not a black hole".into(),
        ));
    }
    Ok(())
}

/// Horizon radii `r± = M ± √(M² − a² − Q²)`.
///
/// # Errors
/// [`RelativityError::InvalidParameter`] for a non-positive mass or a
/// super-extremal hole (no horizon exists).
pub fn horizons(bh: &KerrNewman) -> Result<Horizons> {
    check_blackhole(bh)?;
    let root = bh.horizon_discriminant().max(0.0).sqrt();
    Ok(Horizons {
        outer: bh.mass + root,
        inner: bh.mass - root,
    })
}

/// Outer ergosurface radius at polar angle `theta`:
/// `r_E = M + √(M² − Q² − a² cos²θ)`. Between this surface and the outer
/// horizon lies the ergosphere, where no static observer can exist.
///
/// # Errors
/// [`RelativityError::InvalidParameter`] for an invalid hole.
pub fn ergosphere_radius(bh: &KerrNewman, theta: f64) -> Result<f64> {
    check_blackhole(bh)?;
    let c = theta.cos();
    let disc = bh.mass * bh.mass - bh.charge * bh.charge - bh.spin * bh.spin * c * c;
    Ok(bh.mass + disc.max(0.0).sqrt())
}

/// Radius of the equatorial circular photon orbit (the photon sphere).
///
/// Closed forms:
/// * non-rotating (`a = 0`): `r = (3M + √(9M² − 8Q²)) / 2` — reduces to `3M`
///   for Schwarzschild;
/// * Kerr (`Q = 0`): Bardeen's `r = 2M{1 + cos[⅔ arccos(∓a/M)]}` (upper sign
///   prograde).
///
/// `sense` is consulted only in the rotating case.
///
/// # Errors
/// [`RelativityError::Unsupported`] for the combined spin-and-charge case (no
/// closed form implemented); [`RelativityError::InvalidParameter`] for an
/// invalid hole.
pub fn photon_sphere(bh: &KerrNewman, sense: OrbitSense) -> Result<f64> {
    check_blackhole(bh)?;
    let m = bh.mass;
    if bh.spin == 0.0 {
        let inside = 9.0 * m * m - 8.0 * bh.charge * bh.charge;
        return Ok((3.0 * m + inside.sqrt()) / 2.0);
    }
    if bh.charge == 0.0 {
        let astar = bh.spin / m;
        let arg = -sense.sign() * astar;
        return Ok(2.0 * m * (1.0 + ((2.0 / 3.0) * arg.acos()).cos()));
    }
    Err(RelativityError::Unsupported(
        "photon sphere with both spin and charge: use the geodesic path".into(),
    ))
}

/// Radius of the innermost stable circular orbit (ISCO) for a massive test
/// particle on the equator.
///
/// Closed forms: Schwarzschild `6M`; Kerr via Bardeen–Press–Teukolsky (1972),
/// which gives `6M` at `a = 0` and `M` (prograde) / `9M` (retrograde) at
/// extremal spin.
///
/// # Errors
/// [`RelativityError::Unsupported`] when charge is present (no closed form
/// implemented); [`RelativityError::InvalidParameter`] for an invalid hole.
pub fn isco(bh: &KerrNewman, sense: OrbitSense) -> Result<f64> {
    check_blackhole(bh)?;
    if bh.charge != 0.0 {
        return Err(RelativityError::Unsupported(
            "ISCO with charge: use the geodesic path".into(),
        ));
    }
    let m = bh.mass;
    let astar = bh.spin / m;
    let a2 = astar * astar;
    let z1 = 1.0 + (1.0 - a2).cbrt() * ((1.0 + astar).cbrt() + (1.0 - astar).cbrt());
    let z2 = (3.0 * a2 + z1 * z1).sqrt();
    let r_over_m = 3.0 + z2 - sense.sign() * ((3.0 - z1) * (3.0 + z1 + 2.0 * z2)).sqrt();
    Ok(r_over_m * m)
}

/// Apparent shadow radius (critical impact parameter `b_crit`) of a
/// *non-rotating* black hole, as seen by a distant observer:
/// `b_crit = r_ph / √(1 − 2M/r_ph + Q²/r_ph²)`, where `r_ph` is the photon
/// sphere. For Schwarzschild this is the famous `3√3 M ≈ 5.196 M`.
///
/// # Errors
/// [`RelativityError::Unsupported`] for rotating holes — the Kerr shadow is not
/// a circle and is produced by the ray-tracer; [`RelativityError::InvalidParameter`]
/// for an invalid hole.
pub fn shadow_radius(bh: &KerrNewman) -> Result<f64> {
    check_blackhole(bh)?;
    if bh.spin != 0.0 {
        return Err(RelativityError::Unsupported(
            "closed-form shadow radius only for non-rotating holes; the Kerr \
             shadow is asymmetric and comes from the ray-tracer"
                .into(),
        ));
    }
    let r = photon_sphere(bh, OrbitSense::Prograde)?; // sense irrelevant at a=0
    let (m, q) = (bh.mass, bh.charge);
    let f = 1.0 - 2.0 * m / r + q * q / (r * r);
    Ok(r / f.sqrt())
}

/// Gravitational redshift factor `1 + z = ν_emit / ν_obs` for light emitted by
/// a static observer at equatorial radius `r_emit` and received by one at
/// `r_obs`, using `1 + z = √(−g_tt(r_obs) / −g_tt(r_emit))`. For an observer at
/// infinity this is `1/√(1 − 2M/r_emit)` (Schwarzschild).
///
/// # Errors
/// [`RelativityError::OutsideDomain`] if either radius has no static observer
/// (inside the ergosphere/horizon, where `−g_tt ≤ 0`);
/// [`RelativityError::InvalidParameter`] for an invalid hole.
pub fn gravitational_redshift(bh: &KerrNewman, r_emit: f64, r_obs: f64) -> Result<f64> {
    check_blackhole(bh)?;
    let neg_gtt = |r: f64| -> Result<f64> {
        let g = bh.metric::<f64>([0.0, r, FRAC_PI_2, 0.0]);
        let val = -g[0][0];
        if val <= 0.0 {
            return Err(RelativityError::OutsideDomain(format!(
                "no static observer at r={r} (inside ergosphere/horizon)"
            )));
        }
        Ok(val)
    };
    Ok((neg_gtt(r_obs)? / neg_gtt(r_emit)?).sqrt())
}
