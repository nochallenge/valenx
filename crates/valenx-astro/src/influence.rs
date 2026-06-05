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
        assert!((r - 9.24e8).abs() < 2.0e7, "Earth SOI {r} m (~{:.0} km)", r / 1e3);
    }

    #[test]
    fn earth_hill_sphere_is_about_1_5_million_km() {
        let r = hill_sphere_radius(AU_M, M_EARTH, M_SUN);
        // The textbook value is ≈ 1.5 million km.
        assert!((r - 1.496e9).abs() < 3.0e7, "Earth Hill radius {r} m (~{:.0} km)", r / 1e3);
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
        assert!((by_mass - by_mu).abs() < 1e-3, "only the ratio enters: {by_mass} vs {by_mu}");
    }

    #[test]
    fn non_physical_inputs_return_zero() {
        assert_eq!(sphere_of_influence_radius(-1.0, M_EARTH, M_SUN), 0.0);
        assert_eq!(sphere_of_influence_radius(AU_M, 0.0, M_SUN), 0.0);
        assert_eq!(hill_sphere_radius(AU_M, M_EARTH, f64::NAN), 0.0);
        assert_eq!(hill_sphere_radius(f64::INFINITY, M_EARTH, M_SUN), 0.0);
    }
}
