//! Equations of motion: the per-instant accelerations acting on the
//! vehicle in the planar Earth-centred inertial frame.
//!
//! The three contributions are kept as separate, individually testable
//! functions — inverse-square gravity, aerodynamic drag against the
//! co-rotating atmosphere, and thrust — and summed by [`total_accel`].

use nalgebra::Vector2;

use crate::atmosphere::{self, AtmosphereSample};
use crate::constants::{MU_EARTH, OMEGA_EARTH, R_EARTH};
use crate::vehicle::{DragModel, Vehicle};

/// Velocity of the co-rotating atmosphere at a point (m/s): `ω × r`.
pub fn atmosphere_velocity(position: Vector2<f64>) -> Vector2<f64> {
    Vector2::new(-OMEGA_EARTH * position.y, OMEGA_EARTH * position.x)
}

/// Inverse-square gravitational acceleration (m/s²) toward Earth's
/// centre.
pub fn gravity_accel(position: Vector2<f64>) -> Vector2<f64> {
    let r = position.norm();
    if r < 1.0 {
        return Vector2::zeros();
    }
    -MU_EARTH / (r * r * r) * position
}

/// Geometric altitude above the equatorial radius (m).
pub fn altitude(position: Vector2<f64>) -> f64 {
    position.norm() - R_EARTH
}

/// Sample the atmosphere at the vehicle's current altitude.
pub fn atmosphere_at(position: Vector2<f64>) -> AtmosphereSample {
    atmosphere::sample(altitude(position))
}

/// Aerodynamic drag acceleration (m/s²), opposing the velocity relative
/// to the co-rotating air plus any `wind` (m/s, ECI-frame horizontal).
/// Returns zero in vacuum or at rest.
pub fn drag_accel(
    position: Vector2<f64>,
    velocity: Vector2<f64>,
    mass: f64,
    reference_area: f64,
    drag: &DragModel,
    atmos: &AtmosphereSample,
    wind: Vector2<f64>,
) -> Vector2<f64> {
    if atmos.density <= 0.0 || mass <= 0.0 {
        return Vector2::zeros();
    }
    let v_rel = velocity - atmosphere_velocity(position) - wind;
    let speed = v_rel.norm();
    if speed < 1e-6 {
        return Vector2::zeros();
    }
    let mach = speed / atmos.speed_of_sound;
    let cd = drag.cd(mach);
    let force = 0.5 * atmos.density * speed * speed * cd * reference_area;
    -(force / mass) * (v_rel / speed)
}

/// Dynamic pressure `q = ½ρv_rel²` (Pa), with `v_rel` taken relative to
/// the co-rotating air plus `wind`.
pub fn dynamic_pressure(
    position: Vector2<f64>,
    velocity: Vector2<f64>,
    atmos: &AtmosphereSample,
    wind: Vector2<f64>,
) -> f64 {
    if atmos.density <= 0.0 {
        return 0.0;
    }
    let v_rel = velocity - atmosphere_velocity(position) - wind;
    0.5 * atmos.density * v_rel.norm_squared()
}

/// Thrust acceleration (m/s²) given a thrust magnitude (N) and a unit
/// direction.
pub fn thrust_accel(thrust_newtons: f64, direction: Vector2<f64>, mass: f64) -> Vector2<f64> {
    if mass <= 0.0 || thrust_newtons <= 0.0 {
        return Vector2::zeros();
    }
    (thrust_newtons / mass) * direction
}

/// The full right-hand side: total acceleration from gravity + drag +
/// thrust. `thrust_newtons` already accounts for ambient pressure and
/// throttle, and `direction` is the unit thrust vector from guidance.
#[allow(clippy::too_many_arguments)]
pub fn total_accel(
    position: Vector2<f64>,
    velocity: Vector2<f64>,
    mass: f64,
    vehicle: &Vehicle,
    atmos: &AtmosphereSample,
    thrust_newtons: f64,
    direction: Vector2<f64>,
    wind: Vector2<f64>,
) -> Vector2<f64> {
    gravity_accel(position)
        + drag_accel(
            position,
            velocity,
            mass,
            vehicle.reference_area,
            &vehicle.drag,
            atmos,
            wind,
        )
        + thrust_accel(thrust_newtons, direction, mass)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vehicle::DragModel;

    #[test]
    fn surface_gravity_is_about_9_8() {
        let g = gravity_accel(Vector2::new(R_EARTH, 0.0));
        // Points toward the centre (−x here).
        assert!(g.x < 0.0 && g.y.abs() < 1e-9);
        let mag = g.norm();
        assert!((mag - 9.798).abs() < 0.01, "g = {mag}");
    }

    #[test]
    fn gravity_falls_off_with_square() {
        let g1 = gravity_accel(Vector2::new(R_EARTH, 0.0)).norm();
        let g2 = gravity_accel(Vector2::new(2.0 * R_EARTH, 0.0)).norm();
        // Doubling radius quarters gravity.
        assert!((g1 / g2 - 4.0).abs() < 1e-6);
    }

    #[test]
    fn drag_opposes_relative_velocity_and_vanishes_in_vacuum() {
        let drag = DragModel::generic_launch_vehicle();
        // In atmosphere, moving fast radially outward at low altitude.
        let pos = Vector2::new(R_EARTH + 5_000.0, 0.0);
        let vel = Vector2::new(300.0, 0.0) + atmosphere_velocity(pos);
        let atmos = atmosphere_at(pos);
        let d = drag_accel(pos, vel, 50_000.0, 10.0, &drag, &atmos, Vector2::zeros());
        // Relative velocity is +x, so drag is −x.
        assert!(d.x < 0.0, "drag should oppose motion, got {d:?}");

        // Same speed but above the atmosphere -> no drag.
        let pos_v = Vector2::new(R_EARTH + 200_000.0, 0.0);
        let atmos_v = atmosphere_at(pos_v);
        let dv = drag_accel(
            pos_v,
            vel,
            50_000.0,
            10.0,
            &drag,
            &atmos_v,
            Vector2::zeros(),
        );
        assert_eq!(dv, Vector2::zeros());
    }

    #[test]
    fn wind_increases_dynamic_pressure_and_drag() {
        // At a fixed low-altitude state, a wind opposing the air-relative
        // motion raises both the dynamic pressure and the drag magnitude.
        let drag = DragModel::generic_launch_vehicle();
        let pos = Vector2::new(R_EARTH + 8_000.0, 0.0);
        // Air-relative velocity is +y (downrange) before wind.
        let vel = Vector2::new(0.0, 400.0) + atmosphere_velocity(pos);
        let atmos = atmosphere_at(pos);

        let q0 = dynamic_pressure(pos, vel, &atmos, Vector2::zeros());
        // A head-wind in −y opposes the motion, raising relative speed.
        let head = Vector2::new(0.0, -120.0);
        let q1 = dynamic_pressure(pos, vel, &atmos, head);
        assert!(q1 > q0, "q {q0} -> {q1}");

        let d0 = drag_accel(pos, vel, 50_000.0, 10.0, &drag, &atmos, Vector2::zeros());
        let d1 = drag_accel(pos, vel, 50_000.0, 10.0, &drag, &atmos, head);
        assert!(d1.norm() > d0.norm(), "drag {} -> {}", d0.norm(), d1.norm());
    }

    #[test]
    fn thrust_points_along_direction() {
        let dir = Vector2::new(0.0, 1.0);
        let a = thrust_accel(1_000_000.0, dir, 100_000.0);
        assert!((a.y - 10.0).abs() < 1e-9 && a.x.abs() < 1e-12);
    }

    #[test]
    fn atmosphere_co_rotates_eastward_at_launch() {
        let pos = Vector2::new(R_EARTH, 0.0);
        let v = atmosphere_velocity(pos);
        // Eastward (+y), magnitude ω·R ≈ 465 m/s at the equator.
        assert!(v.x.abs() < 1e-9 && v.y > 0.0);
        assert!((v.y - OMEGA_EARTH * R_EARTH).abs() < 1e-6);
        assert!((v.y - 465.1).abs() < 1.0, "vrot {}", v.y);
    }
}
