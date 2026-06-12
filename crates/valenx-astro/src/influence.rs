//! Gravitational spheres of dominance — the **sphere of influence** and the
//! **Hill sphere** of a body orbiting a heavier primary.
//!
//! Both answer "how far out does this body's gravity hold sway against the
//! primary's?", and both are the standard radii of patched-conic /
//! three-body mission design — but they draw the boundary by different physics:
//!
//! - the **sphere of influence** (Laplace) `r_SOI = a·(m/M)^(2/5)` is where a
//!   probe should *switch* its reference body in a patched-conic trajectory —
//!   inside it the small body's pull is the better primary, outside it the big
//!   one is;
//! - the **Hill sphere** `r_H = a·(m/3M)^(1/3)` is where a satellite can stay
//!   *bound* against the primary's tidal stripping — moons and rings live well
//!   within it (in practice within ~⅓–½ `r_H`).
//!
//! `a` is the orbit radius (semi-major axis) of the small body about the
//! primary, `m` its mass and `M` the primary's. Only the **mass ratio** `m/M`
//! enters, so the two masses may be passed as masses *or* as gravitational
//! parameters `μ = G·M` (the `G` cancels). For a small mass ratio the Hill
//! radius is the larger of the two.

/// Whether the inputs are a physical (positive, finite) orbit radius and mass
/// pair — the radii are otherwise undefined.
fn inputs_ok(orbit_radius_m: f64, body_mass: f64, primary_mass: f64) -> bool {
    orbit_radius_m.is_finite()
        && orbit_radius_m > 0.0
        && body_mass.is_finite()
        && body_mass > 0.0
        && primary_mass.is_finite()
        && primary_mass > 0.0
}

/// The **Laplace sphere of influence** radius (m) of a body of mass `body_mass`
/// orbiting a primary of mass `primary_mass` at orbit radius `orbit_radius_m`:
/// `r_SOI = a·(m/M)^(2/5)`. Inside it a patched-conic trajectory treats the
/// small body as the gravitational primary. `body_mass` and `primary_mass` may
/// be masses or gravitational parameters `μ = G·M` — only their ratio matters.
/// Returns `0` for non-physical inputs (non-finite or non-positive).
pub fn sphere_of_influence_radius(orbit_radius_m: f64, body_mass: f64, primary_mass: f64) -> f64 {
    if !inputs_ok(orbit_radius_m, body_mass, primary_mass) {
        return 0.0;
    }
    orbit_radius_m * (body_mass / primary_mass).powf(0.4)
}

/// The **Hill sphere** radius (m) of a body of mass `body_mass` orbiting a
/// primary of mass `primary_mass` at orbit radius `orbit_radius_m`:
/// `r_H = a·(m/3M)^(1/3)`. A satellite must orbit well within it (in practice
/// `~⅓–½ r_H`) to stay gravitationally bound against the primary's tide.
/// `body_mass` and `primary_mass` may be masses or gravitational parameters
/// `μ = G·M` — only their ratio matters. Returns `0` for non-physical inputs.
pub fn hill_sphere_radius(orbit_radius_m: f64, body_mass: f64, primary_mass: f64) -> f64 {
    if !inputs_ok(orbit_radius_m, body_mass, primary_mass) {
        return 0.0;
    }
    orbit_radius_m * (body_mass / (3.0 * primary_mass)).cbrt()
}

/// The fluid-body Roche coefficient (dimensionless): a fluid (deformable)
/// satellite is tidally disrupted at `d ≈ 2.44·R·(ρ_M/ρ_m)^(1/3)`.
pub const ROCHE_FLUID_COEFFICIENT: f64 = 2.44;

/// Whether the Roche-limit inputs are physical — a positive finite primary
/// radius and a positive finite density pair.
fn roche_inputs_ok(primary_radius_m: f64, primary_density: f64, satellite_density: f64) -> bool {
    primary_radius_m.is_finite()
        && primary_radius_m > 0.0
        && primary_density.is_finite()
        && primary_density > 0.0
        && satellite_density.is_finite()
        && satellite_density > 0.0
}

/// The **rigid-body Roche limit** (m) — the distance from a primary of radius
/// `primary_radius_m` (`R_M`) within which a *rigid* satellite held together only
/// by its own self-gravity is pulled apart by the primary's tide:
/// `d = R_M·(2·ρ_M/ρ_m)^(1/3)`, with `ρ_M = primary_density` and
/// `ρ_m = satellite_density`. Only the density ratio and `R_M` enter (the masses
/// cancel). This is the lower estimate — a real, deformable body disrupts farther
/// out, see [`roche_limit_fluid`]. It complements [`hill_sphere_radius`]: the
/// Hill sphere is where a satellite stays orbitally *bound*, the Roche limit is
/// where its own *body* survives. Returns `0` for non-physical input.
pub fn roche_limit_rigid(
    primary_radius_m: f64,
    primary_density: f64,
    satellite_density: f64,
) -> f64 {
    if !roche_inputs_ok(primary_radius_m, primary_density, satellite_density) {
        return 0.0;
    }
    primary_radius_m * (2.0 * primary_density / satellite_density).cbrt()
}

/// The **fluid-body Roche limit** (m) — the Roche limit for a *fluid*
/// (deformable) satellite, which the tide elongates so it disrupts farther out:
/// `d ≈ 2.44·R_M·(ρ_M/ρ_m)^(1/3)` (see [`ROCHE_FLUID_COEFFICIENT`]). It is the
/// more realistic estimate for rubble-pile moons, and is where planetary ring
/// systems sit (Saturn's rings lie inside it). As with the rigid form only the
/// density ratio and `R_M` enter. Returns `0` for non-physical input.
pub fn roche_limit_fluid(
    primary_radius_m: f64,
    primary_density: f64,
    satellite_density: f64,
) -> f64 {
    if !roche_inputs_ok(primary_radius_m, primary_density, satellite_density) {
        return 0.0;
    }
    ROCHE_FLUID_COEFFICIENT * primary_radius_m * (primary_density / satellite_density).cbrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    // Earth about the Sun — the canonical worked example.
    const AU_M: f64 = 1.495_978_707e11;
    const M_EARTH: f64 = 5.972e24;
    const M_SUN: f64 = 1.989e30;

    #[test]
    fn earth_sphere_of_influence_is_about_924000_km() {
        let r = sphere_of_influence_radius(AU_M, M_EARTH, M_SUN);
        // The textbook value is ≈ 0.924 million km.
        assert!(
            (r - 9.24e8).abs() < 2.0e7,
            "Earth SOI {r} m (~{:.0} km)",
            r / 1e3
        );
    }

    #[test]
    fn earth_hill_sphere_is_about_1_5_million_km() {
        let r = hill_sphere_radius(AU_M, M_EARTH, M_SUN);
        // The textbook value is ≈ 1.5 million km.
        assert!(
            (r - 1.496e9).abs() < 3.0e7,
            "Earth Hill radius {r} m (~{:.0} km)",
            r / 1e3
        );
    }

    #[test]
    fn hill_sphere_exceeds_the_sphere_of_influence_for_a_small_mass_ratio() {
        // For m ≪ M the Hill radius is the larger boundary.
        let soi = sphere_of_influence_radius(AU_M, M_EARTH, M_SUN);
        let hill = hill_sphere_radius(AU_M, M_EARTH, M_SUN);
        assert!(hill > soi, "Hill {hill} should exceed SOI {soi}");
    }

    #[test]
    fn both_radii_scale_linearly_with_orbit_radius() {
        // Both are ∝ a (the mass-ratio factor is fixed), so doubling the orbit
        // radius doubles each radius.
        let s1 = sphere_of_influence_radius(AU_M, M_EARTH, M_SUN);
        let s2 = sphere_of_influence_radius(2.0 * AU_M, M_EARTH, M_SUN);
        assert!((s2 - 2.0 * s1).abs() < 1.0, "SOI ∝ a");
        let h1 = hill_sphere_radius(AU_M, M_EARTH, M_SUN);
        let h2 = hill_sphere_radius(2.0 * AU_M, M_EARTH, M_SUN);
        assert!((h2 - 2.0 * h1).abs() < 1.0, "Hill ∝ a");
    }

    #[test]
    fn only_the_mass_ratio_matters() {
        // Masses or gravitational parameters (μ = G·M) give the same radius,
        // since G cancels in m/M — scale both masses by any factor.
        let g = 6.674_30e-11;
        let by_mass = hill_sphere_radius(AU_M, M_EARTH, M_SUN);
        let by_mu = hill_sphere_radius(AU_M, g * M_EARTH, g * M_SUN);
        assert!(
            (by_mass - by_mu).abs() < 1e-3,
            "only the ratio enters: {by_mass} vs {by_mu}"
        );
    }

    #[test]
    fn non_physical_inputs_return_zero() {
        assert_eq!(sphere_of_influence_radius(-1.0, M_EARTH, M_SUN), 0.0);
        assert_eq!(sphere_of_influence_radius(AU_M, 0.0, M_SUN), 0.0);
        assert_eq!(hill_sphere_radius(AU_M, M_EARTH, f64::NAN), 0.0);
        assert_eq!(hill_sphere_radius(f64::INFINITY, M_EARTH, M_SUN), 0.0);
    }

    #[test]
    fn roche_limit_scales_with_density_ratio_and_primary_radius() {
        // Equal densities: the rigid limit is R·2^(1/3) ≈ 1.26·R, the fluid 2.44·R.
        let r = 1.0_f64;
        let rho = 3000.0_f64;
        let rigid = roche_limit_rigid(r, rho, rho);
        let fluid = roche_limit_fluid(r, rho, rho);
        assert!(
            (rigid - 2.0_f64.cbrt()).abs() < 1e-9,
            "rigid equal-density {rigid}"
        );
        assert!((fluid - 2.44).abs() < 1e-9, "fluid equal-density {fluid}");
        // A fluid (deformable) body is disrupted farther out than a rigid one.
        assert!(
            fluid > rigid,
            "fluid Roche {fluid} should exceed rigid {rigid}"
        );

        // Both scale linearly with the primary radius.
        assert!(
            (roche_limit_rigid(2.0, rho, rho) - 2.0 * rigid).abs() < 1e-9,
            "rigid ∝ R_M"
        );
        assert!(
            (roche_limit_fluid(2.0, rho, rho) - 2.0 * fluid).abs() < 1e-9,
            "fluid ∝ R_M"
        );

        // A denser satellite has a SMALLER Roche limit (∝ ρ_m^(−1/3)): 8× density → ½.
        let denser = roche_limit_rigid(r, rho, 8.0 * rho);
        assert!(
            (denser - rigid / 2.0).abs() < 1e-9,
            "8× satellite density halves it"
        );

        // Reference: the Moon's fluid Roche limit about Earth is ≈ 18,400 km.
        let moon = roche_limit_fluid(6.371e6, 5514.0, 3344.0);
        assert!(
            (moon - 1.84e7).abs() < 6e5,
            "lunar fluid Roche ≈ 18,400 km, got {:.0} km",
            moon / 1e3
        );

        // Non-physical inputs → 0.
        assert_eq!(roche_limit_rigid(-1.0, rho, rho), 0.0);
        assert_eq!(roche_limit_fluid(r, 0.0, rho), 0.0);
        assert_eq!(roche_limit_rigid(r, rho, f64::NAN), 0.0);
    }
}
