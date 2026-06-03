//! Case presets — automotive and aircraft wind-tunnel setups.
//!
//! Setting up a wind-tunnel case correctly takes domain knowledge:
//! a car needs a moving ground, a wing needs the right reference
//! length, the boundary conditions differ. This module bundles that
//! knowledge into named presets so a caller asking for "an automotive
//! run" or "an aircraft cruise run" gets the right setup without
//! re-deriving it.
//!
//! It also provides the **rotating-wheel approximation**
//! ([`rotating_wheel_motion`]): a real spinning wheel imparts a
//! tangential surface velocity to the air. Rather than meshing a
//! rotating boundary, the v1 sets the solid-cell velocity in a wheel
//! region to the local tangential velocity `ω × r` — a direct-forcing
//! approximation that captures the first-order effect (a rotating
//! wheel sheds the flow differently from a stationary one).

use nalgebra::Vector3;

use crate::api::AeroRequest;
use crate::domain::{BoundaryConditions, WindTunnel};
use crate::solver::BodyMotion;
use crate::turbulence::TurbulenceModel;

/// A named wind-tunnel preset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AeroPreset {
    /// A passenger-car run — a moving ground at the road speed,
    /// k-ω SST turbulence, ground-effect-aware sizing.
    AutomotiveRoadCar,
    /// A motorsport / downforce run — like the road car but with a
    /// closer ground and a finer grid to resolve the underbody flow.
    AutomotiveDownforce,
    /// An aircraft cruise run — free-air (no ground), k-ω SST, the
    /// external-aero standard.
    AircraftCruise,
    /// A generic bluff-body run — free-air, k-ε, robust defaults.
    GenericBluffBody,
}

impl AeroPreset {
    /// A human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            AeroPreset::AutomotiveRoadCar => "automotive road car",
            AeroPreset::AutomotiveDownforce => "automotive downforce",
            AeroPreset::AircraftCruise => "aircraft cruise",
            AeroPreset::GenericBluffBody => "generic bluff body",
        }
    }

    /// The boundary-condition set this preset uses.
    ///
    /// `speed` is the free-stream speed — for the automotive presets
    /// the moving ground runs at this speed (the car and the road move
    /// together relative to the still air).
    pub fn boundary(self, speed: f64) -> BoundaryConditions {
        match self {
            AeroPreset::AutomotiveRoadCar | AeroPreset::AutomotiveDownforce => {
                BoundaryConditions::automotive(speed)
            }
            AeroPreset::AircraftCruise | AeroPreset::GenericBluffBody => {
                BoundaryConditions::external_aero()
            }
        }
    }

    /// The turbulence model this preset uses.
    pub fn turbulence(self) -> TurbulenceModel {
        match self {
            AeroPreset::GenericBluffBody => TurbulenceModel::KEpsilon,
            _ => TurbulenceModel::KOmegaSST,
        }
    }

    /// Build a complete [`AeroRequest`] for this preset at the given
    /// free-stream speed.
    pub fn request(self, speed: f64) -> AeroRequest {
        let mut r = AeroRequest::new(speed)
            .with_turbulence(self.turbulence())
            .with_boundary(self.boundary(speed));
        // The automotive presets see road-traffic turbulence;
        // a downforce run wants a finer grid.
        match self {
            AeroPreset::AutomotiveRoadCar => {
                r.turbulence_intensity = 0.02;
            }
            AeroPreset::AutomotiveDownforce => {
                r.turbulence_intensity = 0.02;
                r.sizing.cells_across_body = 24;
            }
            AeroPreset::AircraftCruise => {
                r.turbulence_intensity = 0.005;
            }
            AeroPreset::GenericBluffBody => {}
        }
        r
    }
}

/// Build the solid-cell velocity field for a **rotating wheel**.
///
/// `tunnel` supplies the voxelized body; `axle` is a point on the
/// wheel's axle, `axis` the (unit) axle direction, `omega` the angular
/// speed (rad·s⁻¹), and `region_radius` the radius within which a
/// solid cell is treated as part of the wheel. Every solid cell inside
/// that cylinder gets the tangential velocity `ω·axis × (r − axle)`.
///
/// Returns a [`BodyMotion`] ready to pass to the solver — the wheel's
/// solid cells then spin, while the rest of the body stays static.
pub fn rotating_wheel_motion(
    tunnel: &WindTunnel,
    axle: Vector3<f64>,
    axis: Vector3<f64>,
    omega: f64,
    region_radius: f64,
) -> BodyMotion {
    let g = tunnel.grid;
    let n = g.cell_count();
    let axis = axis.try_normalize(1e-12).unwrap_or_else(Vector3::y);
    let mut solid_velocity = vec![[0.0; 3]; n];

    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                let idx = i + g.nx * (j + g.ny * k);
                if !tunnel.body.is_solid(i, j, k) {
                    continue;
                }
                let (cx, cy, cz) = g.cell_centre(i, j, k);
                let r = Vector3::new(cx, cy, cz) - axle;
                // Distance from the axle line.
                let along = r.dot(&axis);
                let radial = r - along * axis;
                if radial.norm() > region_radius {
                    continue; // outside the wheel region — stays static
                }
                // Tangential velocity v = ω·(axis × radial).
                let v = omega * axis.cross(&radial);
                solid_velocity[idx] = [v.x, v.y, v.z];
            }
        }
    }
    BodyMotion { solid_velocity }
}

/// The angular speed (rad·s⁻¹) of a wheel of radius `wheel_radius`
/// rolling without slipping at road speed `road_speed` — `ω = V/R`.
///
/// A convenience so a caller can give the road speed (which they know)
/// instead of an angular speed (which they would have to compute).
pub fn rolling_wheel_omega(road_speed: f64, wheel_radius: f64) -> f64 {
    if wheel_radius > 1e-9 {
        road_speed / wheel_radius
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::box_body;
    use crate::wind::Wind;

    #[test]
    fn preset_names_are_stable() {
        assert_eq!(AeroPreset::AutomotiveRoadCar.name(), "automotive road car");
        assert_eq!(AeroPreset::AircraftCruise.name(), "aircraft cruise");
    }

    #[test]
    fn automotive_presets_use_a_moving_ground() {
        let bc = AeroPreset::AutomotiveRoadCar.boundary(30.0);
        assert_eq!(bc.z_min, crate::domain::FaceBc::MovingWall(30.0));
        // The aircraft preset is free-air (slip floor).
        let air = AeroPreset::AircraftCruise.boundary(200.0);
        assert_eq!(air.z_min, crate::domain::FaceBc::Slip);
    }

    #[test]
    fn preset_requests_carry_the_right_model() {
        let car = AeroPreset::AutomotiveRoadCar.request(30.0);
        assert_eq!(car.turbulence, TurbulenceModel::KOmegaSST);
        assert!((car.turbulence_intensity - 0.02).abs() < 1e-12);
        let bluff = AeroPreset::GenericBluffBody.request(20.0);
        assert_eq!(bluff.turbulence, TurbulenceModel::KEpsilon);
        // The downforce preset asks for a finer grid.
        let df = AeroPreset::AutomotiveDownforce.request(50.0);
        assert!(df.sizing.cells_across_body >= 24);
    }

    #[test]
    fn rolling_wheel_omega_is_speed_over_radius() {
        // A 0.3 m wheel at 30 m/s → ω = 100 rad/s.
        assert!((rolling_wheel_omega(30.0, 0.3) - 100.0).abs() < 1e-9);
        // Guarded against a zero radius.
        assert_eq!(rolling_wheel_omega(30.0, 0.0), 0.0);
    }

    #[test]
    fn rotating_wheel_motion_spins_the_solid_cells() {
        // A box body, treated as a "wheel" rotating about the y axis:
        // its solid cells must get a non-zero tangential velocity, and
        // the velocity magnitude must scale with the radial distance.
        let body = box_body(
            Vector3::new(-1.0, -0.5, -1.0),
            Vector3::new(1.0, 0.5, 1.0),
        );
        let tunnel = WindTunnel::build(&body, Wind::straight(20.0).unwrap()).unwrap();
        let axle = Vector3::new(0.0, 0.0, 0.0);
        let omega = 50.0;
        let motion = rotating_wheel_motion(
            &tunnel,
            axle,
            Vector3::new(0.0, 1.0, 0.0),
            omega,
            10.0, // large region — covers the whole body
        );
        // Some solid cell must have picked up a tangential velocity.
        let mut max_speed = 0.0f64;
        for v in &motion.solid_velocity {
            max_speed = max_speed.max((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt());
        }
        assert!(max_speed > 0.0, "the wheel should spin its solid cells");
        // The tangential speed at radius r is ω·r — a cell well off
        // the axle moves faster than one near it. The peak speed must
        // be plausible: ω times roughly the body half-diagonal.
        assert!(
            max_speed < omega * 2.0,
            "wheel tangential speed {max_speed} unreasonably high"
        );
    }

    #[test]
    fn rotating_wheel_region_radius_limits_the_spin() {
        // With a tiny region radius, no solid cell qualifies — every
        // cell stays static.
        let body = box_body(
            Vector3::new(-1.0, -1.0, -1.0),
            Vector3::new(1.0, 1.0, 1.0),
        );
        let tunnel = WindTunnel::build(&body, Wind::straight(20.0).unwrap()).unwrap();
        let motion = rotating_wheel_motion(
            &tunnel,
            Vector3::new(100.0, 100.0, 100.0), // axle far away
            Vector3::new(0.0, 1.0, 0.0),
            50.0,
            0.01, // tiny region
        );
        assert!(
            motion
                .solid_velocity
                .iter()
                .all(|v| v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0),
            "no cell should be inside the tiny wheel region"
        );
    }
}
