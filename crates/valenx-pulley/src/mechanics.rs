//! Forces, distances, work and efficiency for a [`PulleySystem`].
//!
//! ## Model
//!
//! For an **ideal** (friction-free, massless) machine the load `W` is
//! shared equally over the `n` supporting rope segments, so the effort is
//!
//! ```text
//! F_ideal = W / MA          (MA = n = ideal mechanical advantage)
//! ```
//!
//! and energy is conserved: the work the operator puts in equals the work
//! done on the load,
//!
//! ```text
//! W_in = F_ideal * s_effort = W * s_load = W_out,
//! ```
//! with the inextensible rope forcing `s_effort = VR * s_load` and
//! `VR = MA`.
//!
//! A **real** machine wastes some input work to sheave-bearing friction
//! and rope stiffness. Efficiency `eta` in `(0, 1]` is the ratio of useful
//! output work to input work, so the operator must apply more force than
//! the ideal:
//!
//! ```text
//! eta   = W_out / W_in
//! F_real = W / (MA * eta) = F_ideal / eta
//! AMA   = W / F_real = MA * eta          (actual mechanical advantage)
//! ```
//!
//! The velocity ratio is a kinematic property of the rope geometry and is
//! unchanged by friction, so `eta = AMA / VR`. The work lost to friction
//! over one lift of the load through `s_load` is
//! `W_in - W_out = W_out (1 / eta - 1)`.
//!
//! ## Honest scope
//!
//! These are idealized rigid-body closed-form relations: the rope is
//! treated as inextensible and massless, sheaves as either friction-free
//! (ideal) or characterised by a single lumped scalar `eta` (real), and
//! the load as a point weight. They reproduce textbook pulley results,
//! NOT the behaviour of real rigging — do not use them to size lifting
//! gear.

use crate::error::{PulleyError, Result};
use crate::system::PulleySystem;

/// Validate that a force / weight magnitude is finite and non-negative.
fn check_load(load: f64) -> Result<()> {
    if !load.is_finite() {
        return Err(PulleyError::invalid("load", "must be a finite number"));
    }
    if load < 0.0 {
        return Err(PulleyError::invalid("load", "must be non-negative"));
    }
    Ok(())
}

/// Validate that an efficiency lies in the half-open interval `(0, 1]`.
fn check_efficiency(eta: f64) -> Result<()> {
    if !eta.is_finite() {
        return Err(PulleyError::invalid(
            "efficiency",
            "must be a finite number",
        ));
    }
    if eta <= 0.0 || eta > 1.0 {
        return Err(PulleyError::invalid(
            "efficiency",
            "must lie in the half-open interval (0, 1]",
        ));
    }
    Ok(())
}

/// The ideal (friction-free) effort needed to hold or slowly raise `load`
/// with `system`:
///
/// ```text
/// F_ideal = load / MA.
/// ```
///
/// More supporting ropes means a larger `MA` and therefore a smaller
/// effort for the same load.
///
/// # Errors
///
/// Returns [`PulleyError::Invalid`] if `load` is negative or non-finite.
pub fn ideal_effort(system: PulleySystem, load: f64) -> Result<f64> {
    check_load(load)?;
    Ok(load / system.ideal_mechanical_advantage())
}

/// The real effort needed once a lumped efficiency `eta` in `(0, 1]` is
/// accounted for:
///
/// ```text
/// F_real = load / (MA * eta) = F_ideal / eta.
/// ```
///
/// Because `eta <= 1`, the real effort is always greater than or equal to
/// the ideal effort, with equality only for the perfect machine
/// `eta = 1`.
///
/// # Errors
///
/// Returns [`PulleyError::Invalid`] if `load` is negative / non-finite or
/// if `eta` is outside `(0, 1]`.
pub fn real_effort(system: PulleySystem, load: f64, eta: f64) -> Result<f64> {
    check_load(load)?;
    check_efficiency(eta)?;
    Ok(load / (system.ideal_mechanical_advantage() * eta))
}

/// The actual mechanical advantage of a real machine,
/// `AMA = load / F_real = MA * eta`. Always less than or equal to the
/// ideal mechanical advantage.
///
/// # Errors
///
/// Returns [`PulleyError::Invalid`] if `eta` is outside `(0, 1]`.
pub fn actual_mechanical_advantage(system: PulleySystem, eta: f64) -> Result<f64> {
    check_efficiency(eta)?;
    Ok(system.ideal_mechanical_advantage() * eta)
}

/// Efficiency recovered from a *measured* actual effort:
///
/// ```text
/// eta = AMA / VR = (load / F_measured) / MA = F_ideal / F_measured.
/// ```
///
/// # Errors
///
/// Returns [`PulleyError::Invalid`] if `load` is negative / non-finite or
/// if `measured_effort` is not strictly positive. Returns
/// [`PulleyError::Inconsistent`] if the measured effort is below the
/// friction-free ideal effort, which would imply `eta > 1`.
pub fn efficiency_from_effort(
    system: PulleySystem,
    load: f64,
    measured_effort: f64,
) -> Result<f64> {
    check_load(load)?;
    if !measured_effort.is_finite() || measured_effort <= 0.0 {
        return Err(PulleyError::invalid(
            "measured_effort",
            "must be a finite positive number",
        ));
    }
    let ideal = load / system.ideal_mechanical_advantage();
    let eta = ideal / measured_effort;
    if eta > 1.0 {
        return Err(PulleyError::inconsistent(
            "measured effort is below the friction-free ideal, implying efficiency > 1",
        ));
    }
    Ok(eta)
}

/// The distance the effort end of the rope must travel to raise the load
/// by `load_distance`:
///
/// ```text
/// s_effort = VR * s_load = MA * s_load.
/// ```
///
/// This is a purely kinematic relation and does not depend on efficiency.
///
/// # Errors
///
/// Returns [`PulleyError::Invalid`] if `load_distance` is negative or
/// non-finite.
pub fn effort_distance(system: PulleySystem, load_distance: f64) -> Result<f64> {
    if !load_distance.is_finite() || load_distance < 0.0 {
        return Err(PulleyError::invalid(
            "load_distance",
            "must be a finite non-negative number",
        ));
    }
    Ok(system.velocity_ratio() * load_distance)
}

/// Useful output work done on the load, `W_out = load * load_distance`.
///
/// # Errors
///
/// Returns [`PulleyError::Invalid`] if `load` or `load_distance` is
/// negative / non-finite.
pub fn output_work(load: f64, load_distance: f64) -> Result<f64> {
    check_load(load)?;
    if !load_distance.is_finite() || load_distance < 0.0 {
        return Err(PulleyError::invalid(
            "load_distance",
            "must be a finite non-negative number",
        ));
    }
    Ok(load * load_distance)
}

/// Input work the operator supplies to raise `load` through
/// `load_distance` with efficiency `eta`:
///
/// ```text
/// W_in = W_out / eta = load * load_distance / eta.
/// ```
///
/// For the ideal machine (`eta = 1`) this equals the output work, so no
/// energy is lost.
///
/// # Errors
///
/// Returns [`PulleyError::Invalid`] if `load` / `load_distance` is
/// negative / non-finite or if `eta` is outside `(0, 1]`.
pub fn input_work(load: f64, load_distance: f64, eta: f64) -> Result<f64> {
    let w_out = output_work(load, load_distance)?;
    check_efficiency(eta)?;
    Ok(w_out / eta)
}

/// Work lost to friction over one lift through `load_distance`:
///
/// ```text
/// W_loss = W_in - W_out = W_out * (1 / eta - 1).
/// ```
///
/// Zero for the ideal machine (`eta = 1`).
///
/// # Errors
///
/// Returns [`PulleyError::Invalid`] if `load` / `load_distance` is
/// negative / non-finite or if `eta` is outside `(0, 1]`.
pub fn work_lost(load: f64, load_distance: f64, eta: f64) -> Result<f64> {
    let w_out = output_work(load, load_distance)?;
    let w_in = input_work(load, load_distance, eta)?;
    Ok(w_in - w_out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// Ideal effort is load / MA, checked against hand-computed values for
    /// each canonical arrangement.
    #[test]
    fn ideal_effort_is_load_over_ma() {
        let load = 600.0;

        // Fixed: MA = 1 -> effort == load.
        let f = ideal_effort(PulleySystem::fixed(), load).unwrap();
        assert!((f - 600.0).abs() < EPS, "got {f}");

        // Movable: MA = 2 -> effort == load / 2.
        let f = ideal_effort(PulleySystem::movable(), load).unwrap();
        assert!((f - 300.0).abs() < EPS, "got {f}");

        // Block-and-tackle, n = 4 -> effort == load / 4.
        let p = PulleySystem::block_and_tackle(4).unwrap();
        let f = ideal_effort(p, load).unwrap();
        assert!((f - 150.0).abs() < EPS, "got {f}");
    }

    /// More supporting ropes means strictly less effort for the same load.
    #[test]
    fn more_ropes_means_less_effort() {
        let load = 1000.0;
        let mut prev = f64::INFINITY;
        for n in 1..=10u32 {
            let p = PulleySystem::block_and_tackle(n).unwrap();
            let f = ideal_effort(p, load).unwrap();
            assert!(
                f < prev,
                "effort did not decrease at n = {n}: {f} !< {prev}"
            );
            prev = f;
        }
    }

    /// Effort * MA == load exactly for the ideal machine (the defining
    /// identity of mechanical advantage).
    #[test]
    fn effort_times_ma_recovers_load() {
        let load = 850.0;
        for n in 1..=12u32 {
            let p = PulleySystem::block_and_tackle(n).unwrap();
            let f = ideal_effort(p, load).unwrap();
            let recovered = f * p.ideal_mechanical_advantage();
            assert!((recovered - load).abs() < EPS, "n = {n}");
        }
    }

    /// Real effort exceeds ideal effort whenever eta < 1, and equals it at
    /// eta == 1.
    #[test]
    fn real_effort_exceeds_ideal_below_unit_efficiency() {
        let p = PulleySystem::block_and_tackle(4).unwrap();
        let load = 400.0;
        let ideal = ideal_effort(p, load).unwrap();

        // eta = 0.8 -> F_real = 100 / 0.8 = 125.
        let real = real_effort(p, load, 0.8).unwrap();
        assert!((real - 125.0).abs() < EPS, "got {real}");
        assert!(real > ideal);

        // eta = 1.0 -> equal.
        let real = real_effort(p, load, 1.0).unwrap();
        assert!((real - ideal).abs() < EPS, "got {real}");
    }

    /// Ideal vs real: actual mechanical advantage is MA * eta, never above
    /// the ideal MA.
    #[test]
    fn actual_ma_is_ideal_times_efficiency() {
        let p = PulleySystem::block_and_tackle(6).unwrap();
        let ama = actual_mechanical_advantage(p, 0.75).unwrap();
        assert!((ama - 4.5).abs() < EPS, "got {ama}");
        assert!(ama <= p.ideal_mechanical_advantage());

        // At eta = 1 the actual MA equals the ideal MA.
        let ama = actual_mechanical_advantage(p, 1.0).unwrap();
        assert!((ama - p.ideal_mechanical_advantage()).abs() < EPS);
    }

    /// Velocity ratio equals MA: effort travels MA times the load travel.
    #[test]
    fn effort_distance_is_ma_times_load_distance() {
        let p = PulleySystem::block_and_tackle(5).unwrap();
        let s = effort_distance(p, 2.0).unwrap();
        assert!((s - 10.0).abs() < EPS, "got {s}");
    }

    /// Ideal machine conserves work: W_in == W_out at eta == 1.
    #[test]
    fn ideal_machine_conserves_work() {
        let p = PulleySystem::block_and_tackle(4).unwrap();
        let load = 200.0;
        let s_load = 3.0;

        // Cross-check the force-distance product against the load-distance
        // product: F_ideal * s_effort == load * s_load.
        let f = ideal_effort(p, load).unwrap();
        let s_effort = effort_distance(p, s_load).unwrap();
        let w_in_force = f * s_effort;
        let w_out = output_work(load, s_load).unwrap();
        assert!((w_in_force - w_out).abs() < EPS, "{w_in_force} != {w_out}");

        // And the explicit input_work at eta = 1 matches.
        let w_in = input_work(load, s_load, 1.0).unwrap();
        assert!((w_in - w_out).abs() < EPS, "{w_in} != {w_out}");
    }

    /// Real machine: input work exceeds output work and the difference is
    /// the friction loss W_out (1/eta - 1).
    #[test]
    fn real_machine_loses_work_to_friction() {
        let load = 500.0;
        let s_load = 2.0;
        let eta = 0.5;

        let w_out = output_work(load, s_load).unwrap(); // 1000 J
        let w_in = input_work(load, s_load, eta).unwrap(); // 2000 J
        let lost = work_lost(load, s_load, eta).unwrap(); // 1000 J

        assert!((w_out - 1000.0).abs() < EPS, "got {w_out}");
        assert!((w_in - 2000.0).abs() < EPS, "got {w_in}");
        assert!((lost - 1000.0).abs() < EPS, "got {lost}");
        assert!((lost - (w_in - w_out)).abs() < EPS);

        // No loss for the ideal machine.
        let lost_ideal = work_lost(load, s_load, 1.0).unwrap();
        assert!(lost_ideal.abs() < EPS, "got {lost_ideal}");
    }

    /// Efficiency recovered from a measured effort is the round-trip
    /// inverse of `real_effort`.
    #[test]
    fn efficiency_round_trips_with_real_effort() {
        let p = PulleySystem::block_and_tackle(4).unwrap();
        let load = 400.0;
        let eta = 0.65;
        let f = real_effort(p, load, eta).unwrap();
        let recovered = efficiency_from_effort(p, load, f).unwrap();
        assert!((recovered - eta).abs() < EPS, "got {recovered}");
    }

    /// efficiency == AMA / VR: a measured effort gives an AMA whose ratio
    /// to the velocity ratio is the efficiency.
    #[test]
    fn efficiency_equals_ama_over_vr() {
        let p = PulleySystem::block_and_tackle(5).unwrap();
        let load = 1000.0;
        let measured = 250.0; // AMA = 1000 / 250 = 4; VR = 5; eta = 0.8.
        let eta = efficiency_from_effort(p, load, measured).unwrap();
        let ama = load / measured;
        let vr = p.velocity_ratio();
        assert!((eta - ama / vr).abs() < EPS, "got {eta}");
        assert!((eta - 0.8).abs() < EPS, "got {eta}");
    }

    /// A measured effort below the friction-free ideal is rejected as
    /// inconsistent (would imply eta > 1).
    #[test]
    fn effort_below_ideal_is_inconsistent() {
        let p = PulleySystem::block_and_tackle(4).unwrap();
        let load = 400.0; // ideal effort = 100.
        let err = efficiency_from_effort(p, load, 80.0).unwrap_err();
        assert_eq!(err.code(), "pulley.inconsistent");
    }

    /// Zero load is admissible and yields zero effort / work everywhere.
    #[test]
    fn zero_load_is_admissible() {
        let p = PulleySystem::movable();
        assert!(ideal_effort(p, 0.0).unwrap().abs() < EPS);
        assert!(real_effort(p, 0.0, 0.5).unwrap().abs() < EPS);
        assert!(output_work(0.0, 5.0).unwrap().abs() < EPS);
        assert!(input_work(0.0, 5.0, 0.5).unwrap().abs() < EPS);
        assert!(work_lost(0.0, 5.0, 0.5).unwrap().abs() < EPS);
    }

    /// Bad inputs are rejected with the right code.
    #[test]
    fn invalid_inputs_rejected() {
        let p = PulleySystem::fixed();

        assert_eq!(ideal_effort(p, -1.0).unwrap_err().code(), "pulley.invalid");
        assert_eq!(
            ideal_effort(p, f64::NAN).unwrap_err().code(),
            "pulley.invalid"
        );

        // Efficiency must be in (0, 1].
        assert_eq!(
            real_effort(p, 10.0, 0.0).unwrap_err().code(),
            "pulley.invalid"
        );
        assert_eq!(
            real_effort(p, 10.0, 1.5).unwrap_err().code(),
            "pulley.invalid"
        );
        assert_eq!(
            real_effort(p, 10.0, -0.2).unwrap_err().code(),
            "pulley.invalid"
        );

        // Measured effort must be strictly positive.
        assert_eq!(
            efficiency_from_effort(p, 10.0, 0.0).unwrap_err().code(),
            "pulley.invalid"
        );

        // Distances must be non-negative.
        assert_eq!(
            effort_distance(p, -1.0).unwrap_err().code(),
            "pulley.invalid"
        );
        assert_eq!(
            output_work(10.0, -1.0).unwrap_err().code(),
            "pulley.invalid"
        );
    }
}
