//! Regenerative-braking energy recovery — how much of a car's kinetic
//! energy a regen-capable drivetrain can put back into the battery, and
//! the driving range that recovered energy buys back.
//!
//! When a car of mass `m` slows from speed `v0` to a lower speed `v1`
//! (both in m/s), the kinetic energy it sheds is the textbook difference
//!
//! ```text
//!   ΔE = ½·m·(v0² − v1²)   (joules)
//! ```
//!
//! A regenerative drivetrain captures only a fraction `η ∈ [0, 1]` of
//! that — the rest is lost to the friction brakes, the motor/inverter,
//! the battery's round-trip inefficiency and rolling/aero drag — so the
//! energy actually returned to the pack is
//!
//! ```text
//!   E_regen = η · ΔE = ½·m·η·(v0² − v1²)
//! ```
//!
//! Because you cannot recover energy while *speeding up*, a request with
//! `v1 > v0` (an acceleration, `ΔE < 0`) yields **zero** recovered energy
//! rather than a non-physical negative "recovery": the friction-brake /
//! regen path is one-directional. The signed energy difference is still
//! available on its own ([`kinetic_energy_delta`]) for callers that want
//! it.
//!
//! Finally, a battery-electric car's *range* is set by how many
//! watt-hours it spends per kilometre, `c` (Wh/km). Converting the
//! recovered joules to watt-hours (`1 Wh = 3600 J`) and dividing by that
//! consumption gives the extra distance the recovered charge restores:
//!
//! ```text
//!   Δrange = E_regen[Wh] / c   (km),   E_regen[Wh] = E_regen[J] / 3600
//! ```
//!
//! All relations here are exact algebra — no integration — so the unit
//! tests pin them straight against these closed forms.
//!
//! ## Honest scope — a preliminary-design energy budget, not a powertrain model
//!
//! This is a **research / preliminary-design grade** lumped-energy
//! accounting tool. It treats regeneration as a single constant
//! round-trip efficiency `η` and range as a single constant consumption
//! `c`. It deliberately does **not** model the motor torque/speed
//! efficiency map, inverter and battery loss as a function of current,
//! state-of-charge or temperature limits on charge acceptance, the blend
//! between regen and friction braking, brake-by-wire dynamics, rolling
//! resistance, aerodynamic drag, grade, or any transient. It makes **no
//! claim of parity** with full vehicle-energy or powertrain tools such as
//! AVL CRUISE, GT-SUITE, Simulink/Simscape, Adams/Car or the like, and is
//! not a substitute for measured drive-cycle (WLTP/EPA) data. Use it for
//! first-order "how much could regen buy back?" estimates only.

use crate::Car;

/// Joules in one watt-hour (`1 Wh = 3600 J`), used to convert recovered
/// braking energy into the watt-hours that drive a range estimate.
pub const JOULES_PER_WH: f64 = 3_600.0;

/// Why a regenerative-braking energy/range computation was rejected.
///
/// Every variant marks an input that would otherwise feed a silent
/// `NaN`/`Inf` into the energy or range algebra (a non-finite or negative
/// mass/speed, an efficiency outside `[0, 1]`, or a non-positive
/// consumption that would divide by zero). Carrying the reason keeps the
/// failure explicit rather than returning a meaningless number.
///
/// This enum is `#[non_exhaustive]`: new variants may be added without it
/// being a breaking change, so downstream `match` arms must include a
/// wildcard.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RegenError {
    /// A mass was non-finite or negative (kg).
    InvalidMass,
    /// A speed was non-finite or negative (m/s).
    InvalidSpeed,
    /// The recovery efficiency was non-finite or outside the closed
    /// interval `[0, 1]`.
    InvalidEfficiency,
    /// The energy consumption was non-finite or not strictly positive
    /// (Wh/km), which would make the range division a `NaN`/`Inf`.
    InvalidConsumption,
}

impl core::fmt::Display for RegenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let msg = match self {
            RegenError::InvalidMass => "mass must be finite and >= 0 (kg)",
            RegenError::InvalidSpeed => "speed must be finite and >= 0 (m/s)",
            RegenError::InvalidEfficiency => "efficiency must be finite and in [0, 1]",
            RegenError::InvalidConsumption => "consumption must be finite and > 0 (Wh/km)",
        };
        f.write_str(msg)
    }
}

impl std::error::Error for RegenError {}

/// Shorthand for a regen result that may fail validation.
pub type Result<T> = core::result::Result<T, RegenError>;

/// Validate a mass (kg): finite and non-negative.
fn check_mass(mass: f64) -> Result<()> {
    if !mass.is_finite() || mass < 0.0 {
        return Err(RegenError::InvalidMass);
    }
    Ok(())
}

/// Validate a speed (m/s): finite and non-negative.
fn check_speed(v: f64) -> Result<()> {
    if !v.is_finite() || v < 0.0 {
        return Err(RegenError::InvalidSpeed);
    }
    Ok(())
}

/// Validate an efficiency: finite and within the closed interval `[0, 1]`.
fn check_efficiency(eta: f64) -> Result<()> {
    if !eta.is_finite() || !(0.0..=1.0).contains(&eta) {
        return Err(RegenError::InvalidEfficiency);
    }
    Ok(())
}

/// The **signed** change in kinetic energy (joules) when a body of mass
/// `mass` (kg) changes speed from `v0` to `v1` (both m/s):
///
/// ```text
///   ΔE = ½·m·(v0² − v1²)
/// ```
///
/// This is positive when the body **slows** (`v0 > v1`, energy released),
/// zero at constant speed, and **negative** when it **speeds up**
/// (`v1 > v0`, energy must be supplied). For the one-directional energy
/// that regen can actually recover use [`recoverable_braking_energy`].
///
/// # Errors
///
/// Returns [`RegenError::InvalidMass`] if `mass` is non-finite or
/// negative, or [`RegenError::InvalidSpeed`] if either speed is
/// non-finite or negative.
pub fn kinetic_energy_delta(mass: f64, v0: f64, v1: f64) -> Result<f64> {
    check_mass(mass)?;
    check_speed(v0)?;
    check_speed(v1)?;
    Ok(0.5 * mass * (v0 * v0 - v1 * v1))
}

/// The kinetic energy (joules) **recoverable** by a regenerative
/// drivetrain of round-trip efficiency `efficiency` when a body of mass
/// `mass` (kg) decelerates from `v0` to `v1` (m/s):
///
/// ```text
///   E_regen = η · ½·m·(v0² − v1²),   η = efficiency ∈ [0, 1]
/// ```
///
/// Regeneration is one-directional: if the request is actually an
/// acceleration (`v1 > v0`, so the kinetic-energy change is negative)
/// this returns **0.0**, never a negative value — you cannot harvest
/// energy while speeding up. Consequently the result is always in the
/// closed interval `[0, ½·m·(v0² − v1²)]`: it is `0` at `η = 0`, the full
/// shed kinetic energy at `η = 1`, and otherwise strictly less than it.
///
/// # Errors
///
/// Returns [`RegenError::InvalidMass`] / [`RegenError::InvalidSpeed`] for
/// a non-physical mass or speed, or [`RegenError::InvalidEfficiency`] if
/// `efficiency` is non-finite or outside `[0, 1]`.
pub fn recoverable_braking_energy(mass: f64, v0: f64, v1: f64, efficiency: f64) -> Result<f64> {
    check_efficiency(efficiency)?;
    let delta = kinetic_energy_delta(mass, v0, v1)?;
    // One-directional: clamp an acceleration (ΔE < 0) to zero recovery.
    let shed = delta.max(0.0);
    Ok(efficiency * shed)
}

/// Convert an energy in joules to watt-hours (`Wh = J / 3600`).
///
/// # Errors
///
/// Returns [`RegenError::InvalidConsumption`] only as a stand-in for a
/// non-finite input; a non-finite `joules` is rejected so the conversion
/// cannot silently propagate a `NaN`/`Inf` into a range estimate. (Energy
/// itself is otherwise unconstrained in sign here.)
pub fn joules_to_wh(joules: f64) -> Result<f64> {
    if !joules.is_finite() {
        return Err(RegenError::InvalidConsumption);
    }
    Ok(joules / JOULES_PER_WH)
}

/// The extra driving range (km) that regeneratively recovering the
/// braking energy of a `v0 → v1` deceleration adds back to a
/// battery-electric car, given its energy consumption `consumption_wh_per_km`
/// (Wh/km):
///
/// ```text
///   Δrange = E_regen[Wh] / c,
///   E_regen[Wh] = η·½·m·(v0² − v1²) / 3600
/// ```
///
/// As with [`recoverable_braking_energy`], an acceleration request
/// (`v1 > v0`) recovers nothing and so adds `0.0` km of range.
///
/// # Errors
///
/// Returns the mass / speed / efficiency errors of
/// [`recoverable_braking_energy`], or [`RegenError::InvalidConsumption`]
/// if `consumption_wh_per_km` is non-finite or not strictly positive
/// (which would make the division a `NaN`/`Inf`).
pub fn regen_range_added_km(
    mass: f64,
    v0: f64,
    v1: f64,
    efficiency: f64,
    consumption_wh_per_km: f64,
) -> Result<f64> {
    if !consumption_wh_per_km.is_finite() || consumption_wh_per_km <= 0.0 {
        return Err(RegenError::InvalidConsumption);
    }
    let e_j = recoverable_braking_energy(mass, v0, v1, efficiency)?;
    let e_wh = e_j / JOULES_PER_WH;
    Ok(e_wh / consumption_wh_per_km)
}

/// A small regenerative-braking model: a constant round-trip recovery
/// efficiency plus a constant battery energy consumption, applied to a
/// [`Car`] by reading its mass.
///
/// Construct it with [`RegenBraking::new`] (which validates the
/// parameters once) and then ask it, for any deceleration of a given
/// `Car`, how much energy comes back ([`RegenBraking::recovered_energy_j`])
/// and how much range that buys ([`RegenBraking::range_added_km`]). The
/// car's kinetic energy is taken from its `mass` field; this type holds
/// only the drivetrain-side parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RegenBraking {
    /// Round-trip regeneration efficiency `η ∈ [0, 1]` (dimensionless).
    pub efficiency: f64,
    /// Battery energy consumption used for the range estimate (Wh/km).
    pub consumption_wh_per_km: f64,
}

impl RegenBraking {
    /// Build a validated regen model from a recovery `efficiency`
    /// (`∈ [0, 1]`) and a battery `consumption_wh_per_km` (`> 0`).
    ///
    /// # Errors
    ///
    /// Returns [`RegenError::InvalidEfficiency`] if `efficiency` is
    /// non-finite or outside `[0, 1]`, or
    /// [`RegenError::InvalidConsumption`] if `consumption_wh_per_km` is
    /// non-finite or not strictly positive.
    pub fn new(efficiency: f64, consumption_wh_per_km: f64) -> Result<Self> {
        check_efficiency(efficiency)?;
        if !consumption_wh_per_km.is_finite() || consumption_wh_per_km <= 0.0 {
            return Err(RegenError::InvalidConsumption);
        }
        Ok(Self {
            efficiency,
            consumption_wh_per_km,
        })
    }

    /// Energy (joules) this drivetrain recovers when `car` decelerates
    /// from `v0` to `v1` (m/s). Reads the car's mass; see
    /// [`recoverable_braking_energy`] for the closed form and the
    /// zero-on-acceleration guard.
    ///
    /// # Errors
    ///
    /// Returns [`RegenError::InvalidMass`] if the car's mass is
    /// non-finite or negative, or [`RegenError::InvalidSpeed`] for a
    /// non-physical speed.
    pub fn recovered_energy_j(&self, car: &Car, v0: f64, v1: f64) -> Result<f64> {
        recoverable_braking_energy(car.mass, v0, v1, self.efficiency)
    }

    /// Extra driving range (km) this drivetrain adds back by recovering
    /// the braking energy of a `v0 → v1` deceleration of `car`. Reads the
    /// car's mass; see [`regen_range_added_km`] for the closed form.
    ///
    /// # Errors
    ///
    /// Returns [`RegenError::InvalidMass`] / [`RegenError::InvalidSpeed`]
    /// for a non-physical car mass or speed.
    pub fn range_added_km(&self, car: &Car, v0: f64, v1: f64) -> Result<f64> {
        regen_range_added_km(
            car.mass,
            v0,
            v1,
            self.efficiency,
            self.consumption_wh_per_km,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A representative compact-EV mass for the worked numbers below.
    const M: f64 = 1_500.0; // kg

    #[test]
    fn full_stop_recovers_half_m_v0_sq_times_eff() {
        // ORACLE: stopping completely (v1 = 0) sheds the entire kinetic
        // energy ½·m·v0², of which η is recovered.
        // m = 1500 kg, v0 = 30 m/s, η = 0.6:
        //   ½·1500·30² = 0.5·1500·900 = 675_000 J
        //   E_regen   = 0.6·675_000  = 405_000 J
        let v0 = 30.0;
        let eta = 0.6;
        let e = recoverable_braking_energy(M, v0, 0.0, eta).expect("valid");
        let expected = 0.5 * M * v0 * v0 * eta;
        assert!((e - expected).abs() < 1e-9, "E_regen = {e}");
        assert!((e - 405_000.0).abs() < 1e-6, "E_regen = {e}");
    }

    #[test]
    fn recovered_never_exceeds_kinetic_shed() {
        // E_regen = η·ΔE with η ∈ [0,1] and ΔE >= 0, so it can never
        // exceed the kinetic energy released by the deceleration.
        let v0 = 40.0;
        let v1 = 12.0;
        let delta = kinetic_energy_delta(M, v0, v1).expect("valid");
        assert!(delta > 0.0);
        for &eta in &[0.0, 0.25, 0.5, 0.85, 1.0] {
            let e = recoverable_braking_energy(M, v0, v1, eta).expect("valid");
            assert!(e <= delta + 1e-9, "eta={eta}: {e} > shed {delta}");
            assert!(e >= -1e-12, "eta={eta}: negative recovery {e}");
        }
    }

    #[test]
    fn efficiency_zero_recovers_nothing() {
        // η = 0 ⇒ E_regen = 0 regardless of the speeds.
        let e = recoverable_braking_energy(M, 33.0, 5.0, 0.0).expect("valid");
        assert_eq!(e, 0.0, "eff=0 must recover exactly nothing");
    }

    #[test]
    fn efficiency_one_recovers_full_kinetic_delta() {
        // η = 1 ⇒ E_regen equals the full kinetic-energy difference
        // ½·m·(v0² − v1²) exactly.
        let v0 = 33.0;
        let v1 = 5.0;
        let e = recoverable_braking_energy(M, v0, v1, 1.0).expect("valid");
        let delta = kinetic_energy_delta(M, v0, v1).expect("valid");
        let expected = 0.5 * M * (v0 * v0 - v1 * v1);
        assert!((e - delta).abs() < 1e-9, "eff=1 vs ΔE: {e} vs {delta}");
        assert!((e - expected).abs() < 1e-9, "eff=1 closed form: {e}");
    }

    #[test]
    fn acceleration_request_recovers_zero() {
        // v1 > v0 is an acceleration: ΔE is negative, but recoverable
        // energy is guarded to exactly 0 (you cannot regen while
        // speeding up). The signed delta itself stays negative.
        let v0 = 10.0;
        let v1 = 25.0;
        let delta = kinetic_energy_delta(M, v0, v1).expect("valid");
        assert!(delta < 0.0, "signed ΔE should be negative, got {delta}");
        // Even at full efficiency, recovery clamps to zero — not negative.
        let e = recoverable_braking_energy(M, v0, v1, 1.0).expect("valid");
        assert_eq!(e, 0.0, "accel must recover zero, got {e}");
        // ...and so does the range it would add.
        let dr = regen_range_added_km(M, v0, v1, 1.0, 150.0).expect("valid");
        assert_eq!(dr, 0.0, "accel must add zero range, got {dr}");
    }

    #[test]
    fn range_added_matches_closed_form() {
        // ORACLE: full stop from v0 = 30 m/s, m = 1500 kg, η = 0.6,
        // consumption c = 150 Wh/km.
        //   E_regen = 405_000 J  (see full_stop test)
        //   in Wh    = 405_000 / 3600 = 112.5 Wh
        //   Δrange   = 112.5 / 150   = 0.75 km
        let dr = regen_range_added_km(M, 30.0, 0.0, 0.6, 150.0).expect("valid");
        let e_wh = (0.5 * M * 30.0 * 30.0 * 0.6) / JOULES_PER_WH;
        let expected = e_wh / 150.0;
        assert!((dr - expected).abs() < 1e-12, "Δrange = {dr}");
        assert!((dr - 0.75).abs() < 1e-9, "Δrange = {dr} km");
    }

    #[test]
    fn joules_to_wh_is_exact() {
        // 3_600_000 J = 1000 Wh = 1 kWh.
        let wh = joules_to_wh(3_600_000.0).expect("valid");
        assert!((wh - 1_000.0).abs() < 1e-9, "Wh = {wh}");
    }

    #[test]
    fn rejects_non_physical_inputs() {
        // Mass.
        assert_eq!(
            recoverable_braking_energy(-1.0, 10.0, 0.0, 0.5),
            Err(RegenError::InvalidMass)
        );
        assert_eq!(
            recoverable_braking_energy(f64::NAN, 10.0, 0.0, 0.5),
            Err(RegenError::InvalidMass)
        );
        // Speed.
        assert_eq!(
            recoverable_braking_energy(M, -1.0, 0.0, 0.5),
            Err(RegenError::InvalidSpeed)
        );
        assert_eq!(
            recoverable_braking_energy(M, 10.0, f64::INFINITY, 0.5),
            Err(RegenError::InvalidSpeed)
        );
        // Efficiency outside [0, 1].
        assert_eq!(
            recoverable_braking_energy(M, 10.0, 0.0, 1.5),
            Err(RegenError::InvalidEfficiency)
        );
        assert_eq!(
            recoverable_braking_energy(M, 10.0, 0.0, -0.1),
            Err(RegenError::InvalidEfficiency)
        );
        // Consumption.
        assert_eq!(
            regen_range_added_km(M, 30.0, 0.0, 0.6, 0.0),
            Err(RegenError::InvalidConsumption)
        );
        assert_eq!(
            regen_range_added_km(M, 30.0, 0.0, 0.6, -5.0),
            Err(RegenError::InvalidConsumption)
        );
    }

    #[test]
    fn regen_braking_new_validates() {
        assert!(RegenBraking::new(0.6, 150.0).is_ok());
        assert_eq!(
            RegenBraking::new(1.2, 150.0),
            Err(RegenError::InvalidEfficiency)
        );
        assert_eq!(
            RegenBraking::new(0.6, 0.0),
            Err(RegenError::InvalidConsumption)
        );
        assert_eq!(
            RegenBraking::new(0.6, f64::NAN),
            Err(RegenError::InvalidConsumption)
        );
    }
}
