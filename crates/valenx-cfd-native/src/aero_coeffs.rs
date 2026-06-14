//! **Dimensionless aerodynamic coefficient helpers** — the textbook
//! scalar conversions between a measured force and the non-dimensional
//! numbers that report it: dynamic pressure, the drag and lift
//! coefficients, and the Reynolds number.
//!
//! # What's here
//!
//! Four pure functions, each the standard closed-form definition. Given
//! a free-stream density `ρ` and speed `U`:
//!
//! - **Dynamic pressure** `q = ½·ρ·U²` ([`dynamic_pressure`]) — the
//!   kinetic energy per unit volume of the free stream, the pressure
//!   scale every aerodynamic coefficient normalises by.
//! - **Drag coefficient** `C_d = 2·F_d / (ρ·U²·A) = F_d / (q·A)`
//!   ([`drag_coefficient`]) — the streamwise force `F_d` made
//!   dimensionless by the dynamic pressure and a reference area `A`.
//! - **Lift coefficient** `C_l = 2·F_l / (ρ·U²·A) = F_l / (q·A)`
//!   ([`lift_coefficient`]) — the cross-stream force `F_l`, normalised
//!   identically. (`C_d` and `C_l` differ only in which force component
//!   you feed them, so the two share one definition.)
//! - **Reynolds number** `Re = ρ·U·L / μ` ([`reynolds`]) — the ratio of
//!   inertial to viscous forces over a reference length `L`, with `μ`
//!   the *dynamic* viscosity. (For the *kinematic* viscosity `ν = μ/ρ`,
//!   `Re = U·L / ν`.)
//!
//! These are the bookkeeping that turns a solver's (or a wind-tunnel's)
//! raw force in newtons into the `C_d` / `C_l` / `Re` triple a
//! preliminary-design study reports, and back again.
//!
//! # Divide-by-zero policy
//!
//! Each function guards its denominator: a zero (or non-finite) density,
//! speed, area, length, or viscosity that would otherwise produce a
//! `NaN` or an infinity yields **`0.0`** instead. This keeps a coefficient
//! sweep that brackets `U = 0` (a stagnation / start-up point) finite and
//! plottable rather than seeding `NaN` through the rest of a computation.
//! The guard is a deliberate, documented convention — not a physical
//! claim that the coefficient *is* zero there.
//!
//! # Honest scope — definitions, not a force model
//!
//! This module is **research / preliminary-design-grade dimensional
//! bookkeeping**, not an aerodynamics package. It only rescales a force
//! you already have; it computes **no** forces itself — there is no
//! panel method, no boundary-layer / skin-friction model, no
//! compressibility (Mach / Prandtl-Glauert) correction, and no notion of
//! what a "reference area" should be for a given body (planform, frontal,
//! and wetted areas all differ and are the caller's choice). The numbers
//! are exactly as good as the force, density, speed, length, area, and
//! viscosity fed in, all in one consistent unit system (SI throughout
//! here). For validated force prediction use a real CFD run (the native
//! SIMPLE solver in this crate, or Valenx's OpenFOAM / SU2 adapters); do
//! not mistake these conversions for parity with Ansys Fluent or a
//! certified aerodynamic database.

/// Free-stream **dynamic pressure** `q = ½·ρ·U²` (Pa, for SI inputs).
///
/// The kinetic energy per unit volume of a fluid of density `rho`
/// (kg/m³) moving at speed `u` (m/s) — the pressure scale the drag and
/// lift coefficients normalise by.
///
/// Sign note: `u` enters squared, so the result is independent of flow
/// direction and is `≥ 0` for any real input.
///
/// # Examples
///
/// ```
/// use valenx_cfd_native::aero_coeffs::dynamic_pressure;
/// // Sea-level air (ρ ≈ 1.225 kg/m³) at 100 m/s.
/// let q = dynamic_pressure(1.225, 100.0);
/// assert!((q - 6125.0).abs() < 1e-9);
/// ```
#[must_use]
pub fn dynamic_pressure(rho: f64, u: f64) -> f64 {
    let q = 0.5 * rho * u * u;
    if q.is_finite() {
        q
    } else {
        0.0
    }
}

/// **Drag coefficient** `C_d = 2·F_d / (ρ·U²·A)` from a measured
/// streamwise force.
///
/// `force` is the drag force `F_d` (N), `rho` the free-stream density
/// (kg/m³), `u` the free-stream speed (m/s), and `area` the reference
/// area `A` (m²). Equivalently `C_d = F_d / (q·A)` with the dynamic
/// pressure `q = ½·ρ·U²` from [`dynamic_pressure`].
///
/// Returns `0.0` when the denominator `ρ·U²·A` is zero or non-finite
/// (e.g. `u = 0` or `area = 0`) rather than a `NaN` / infinity — see the
/// module-level divide-by-zero policy.
///
/// # Examples
///
/// ```
/// use valenx_cfd_native::aero_coeffs::{drag_coefficient, dynamic_pressure};
/// // Recover a known C_d from the force it implies.
/// let (rho, u, area, cd) = (1.225, 50.0, 2.0, 0.30);
/// let force = cd * dynamic_pressure(rho, u) * area;
/// let back = drag_coefficient(force, rho, u, area);
/// assert!((back - cd).abs() < 1e-12);
/// ```
#[must_use]
pub fn drag_coefficient(force: f64, rho: f64, u: f64, area: f64) -> f64 {
    let denom = rho * u * u * area;
    if denom == 0.0 || !denom.is_finite() {
        return 0.0;
    }
    let cd = 2.0 * force / denom;
    if cd.is_finite() {
        cd
    } else {
        0.0
    }
}

/// **Lift coefficient** `C_l = 2·F_l / (ρ·U²·A)` from a measured
/// cross-stream force.
///
/// Identical in form to [`drag_coefficient`] — only the force component
/// differs: `force` here is the lift force `F_l` (N), normal to the free
/// stream. `rho` is density (kg/m³), `u` the free-stream speed (m/s), and
/// `area` the reference area `A` (m²). Equivalently `C_l = F_l / (q·A)`.
///
/// Returns `0.0` when the denominator `ρ·U²·A` is zero or non-finite,
/// per the module-level divide-by-zero policy.
///
/// # Examples
///
/// ```
/// use valenx_cfd_native::aero_coeffs::{lift_coefficient, dynamic_pressure};
/// let (rho, u, area, cl) = (1.225, 60.0, 1.5, 0.80);
/// let force = cl * dynamic_pressure(rho, u) * area;
/// assert!((lift_coefficient(force, rho, u, area) - cl).abs() < 1e-12);
/// ```
#[must_use]
pub fn lift_coefficient(force: f64, rho: f64, u: f64, area: f64) -> f64 {
    // Lift and drag coefficients share the identical normalisation; they
    // differ only in which force component the caller supplies.
    drag_coefficient(force, rho, u, area)
}

/// **Reynolds number** `Re = ρ·U·L / μ` over a reference length.
///
/// `rho` is the fluid density (kg/m³), `u` the characteristic speed
/// (m/s), `length` the reference length `L` (m), and `mu` the **dynamic**
/// viscosity (Pa·s). The result is dimensionless — the ratio of inertial
/// to viscous forces.
///
/// For the **kinematic** viscosity `ν = μ/ρ` (m²/s) the same number is
/// `Re = U·L / ν`; pass `rho = 1.0` and `mu = ν` to evaluate that form,
/// matching the `ν = 1/Re` convention the rest of this crate's
/// benchmarks use.
///
/// Returns `0.0` when the denominator `μ` is zero or non-finite, per the
/// module-level divide-by-zero policy.
///
/// # Examples
///
/// ```
/// use valenx_cfd_native::aero_coeffs::reynolds;
/// // ρ=1, U=2, L=3, μ=0.5  →  Re = 1·2·3 / 0.5 = 12.
/// assert!((reynolds(1.0, 2.0, 3.0, 0.5) - 12.0).abs() < 1e-12);
/// ```
#[must_use]
pub fn reynolds(rho: f64, u: f64, length: f64, mu: f64) -> f64 {
    if mu == 0.0 || !mu.is_finite() {
        return 0.0;
    }
    let re = rho * u * length / mu;
    if re.is_finite() {
        re
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dynamic_pressure_is_half_rho_u_squared() {
        // Closed form q = ½·ρ·U². ρ=1.225, U=100 → 0.5·1.225·10000 = 6125.
        let q = dynamic_pressure(1.225, 100.0);
        assert!(
            (q - 6125.0).abs() < 1e-9,
            "q {q} should equal 0.5*1.225*100^2 = 6125"
        );
        // Direction-independent: U and -U give the same q.
        assert!((dynamic_pressure(1.225, -100.0) - q).abs() < 1e-12);
        // A second independent point: ρ=2, U=3 → 0.5·2·9 = 9.
        assert!((dynamic_pressure(2.0, 3.0) - 9.0).abs() < 1e-12);
    }

    #[test]
    fn drag_coefficient_inverts_the_drag_formula() {
        // The round-trip identity F = C_d · q · A must recover C_d
        // exactly, where q = ½·ρ·U². This is the defining inverse
        // relationship of the coefficient.
        let rho = 1.225;
        let u = 50.0;
        let area = 2.0;
        let cd_in = 0.30;
        let q = dynamic_pressure(rho, u);
        let force = cd_in * q * area; // F = C_d · q · A
        let cd_out = drag_coefficient(force, rho, u, area);
        assert!(
            (cd_out - cd_in).abs() < 1e-12,
            "C_d round-trip: got {cd_out}, expected {cd_in}"
        );
        // And the explicit 2F/(ρU²A) form agrees with F/(qA).
        let direct = 2.0 * force / (rho * u * u * area);
        assert!((cd_out - direct).abs() < 1e-12);
    }

    #[test]
    fn lift_coefficient_inverts_the_lift_formula() {
        // Same normalisation as drag, fed the cross-stream force.
        let rho = 1.225;
        let u = 60.0;
        let area = 1.5;
        let cl_in = 0.80;
        let force = cl_in * dynamic_pressure(rho, u) * area;
        let cl_out = lift_coefficient(force, rho, u, area);
        assert!(
            (cl_out - cl_in).abs() < 1e-12,
            "C_l round-trip: got {cl_out}, expected {cl_in}"
        );
        // Lift and drag share the formula: same args → same number.
        assert!(
            (lift_coefficient(force, rho, u, area) - drag_coefficient(force, rho, u, area)).abs()
                < 1e-15
        );
    }

    #[test]
    fn reynolds_matches_known_case() {
        // ρ=1, U=2, L=3, μ=0.5 → Re = 1·2·3 / 0.5 = 12 (closed form).
        let re = reynolds(1.0, 2.0, 3.0, 0.5);
        assert!((re - 12.0).abs() < 1e-12, "Re {re} should equal 12");

        // A physical case: 20 °C water-like flow, ρ=1000, U=1 m/s over
        // L=0.1 m, μ=1.0e-3 Pa·s → Re = 1000·1·0.1 / 1e-3 = 1.0e5.
        let re_water = reynolds(1000.0, 1.0, 0.1, 1.0e-3);
        assert!(
            (re_water - 1.0e5).abs() < 1e-6,
            "water Re {re_water} should equal 1.0e5"
        );

        // Kinematic form Re = U·L/ν via rho=1, mu=ν: U=2, L=3, ν=0.5 → 12.
        assert!((reynolds(1.0, 2.0, 3.0, 0.5) - 12.0).abs() < 1e-12);
    }

    #[test]
    fn degenerate_inputs_return_zero_not_nan() {
        // Zero speed → coefficient denominators vanish → 0.0, finite.
        let cd_u0 = drag_coefficient(10.0, 1.225, 0.0, 2.0);
        assert_eq!(cd_u0, 0.0);
        assert!(cd_u0.is_finite());

        // Zero reference area → 0.0, finite.
        let cl_a0 = lift_coefficient(10.0, 1.225, 50.0, 0.0);
        assert_eq!(cl_a0, 0.0);
        assert!(cl_a0.is_finite());

        // Zero viscosity → Reynolds guard returns 0.0, not infinity/NaN.
        let re_mu0 = reynolds(1.0, 2.0, 3.0, 0.0);
        assert_eq!(re_mu0, 0.0);
        assert!(re_mu0.is_finite());

        // Zero density (vacuum) also collapses the coefficient denom.
        let cd_rho0 = drag_coefficient(10.0, 0.0, 50.0, 2.0);
        assert_eq!(cd_rho0, 0.0);

        // Dynamic pressure at U=0 is genuinely 0 (and finite).
        let q0 = dynamic_pressure(1.225, 0.0);
        assert_eq!(q0, 0.0);
        assert!(q0.is_finite());

        // Non-finite guards: a NaN viscosity must not poison the result.
        assert_eq!(reynolds(1.0, 2.0, 3.0, f64::NAN), 0.0);
        assert_eq!(drag_coefficient(10.0, 1.225, f64::INFINITY, 2.0), 0.0);
    }
}
