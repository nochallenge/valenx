//! Wind specification and air properties — the free-stream the
//! virtual wind tunnel blows over the body.
//!
//! The [`Wind`] struct bundles everything that defines the on-coming
//! flow: the free-stream speed, its direction (set as a yaw / pitch
//! pair so a caller thinks in aircraft / vehicle angles, not raw
//! vectors), the air density and viscosity, and the upstream
//! turbulence intensity. [`Air`] carries the constant fluid
//! properties; the [`Air::sea_level`] constructor is standard-
//! atmosphere air at 15 °C.

use nalgebra::Vector3;

use crate::error::AeroError;

/// Constant air properties.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Air {
    /// Density `ρ` (kg·m⁻³).
    pub density: f64,
    /// Dynamic viscosity `μ` (Pa·s).
    pub dynamic_viscosity: f64,
}

impl Air {
    /// Build an air model; both properties must be positive.
    pub fn new(density: f64, dynamic_viscosity: f64) -> Result<Air, AeroError> {
        if !(density.is_finite() && density > 0.0) {
            return Err(AeroError::BadParameter {
                name: "density",
                reason: format!("must be positive and finite, got {density}"),
            });
        }
        if !(dynamic_viscosity.is_finite() && dynamic_viscosity > 0.0) {
            return Err(AeroError::BadParameter {
                name: "dynamic_viscosity",
                reason: format!("must be positive and finite, got {dynamic_viscosity}"),
            });
        }
        Ok(Air {
            density,
            dynamic_viscosity,
        })
    }

    /// Standard-atmosphere sea-level air at 15 °C:
    /// `ρ = 1.225 kg·m⁻³`, `μ = 1.81e-5 Pa·s`.
    pub fn sea_level() -> Air {
        Air {
            density: 1.225,
            dynamic_viscosity: 1.81e-5,
        }
    }

    /// The kinematic viscosity `ν = μ / ρ` (m²·s⁻¹).
    #[inline]
    pub fn kinematic_viscosity(&self) -> f64 {
        self.dynamic_viscosity / self.density
    }
}

impl Default for Air {
    /// Sea-level standard air.
    fn default() -> Self {
        Air::sea_level()
    }
}

/// The on-coming free-stream specification for a wind-tunnel run.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Wind {
    /// Free-stream speed `U∞` (m·s⁻¹).
    pub speed: f64,
    /// Yaw angle (radians) — rotation of the wind about the vertical
    /// `z` axis. Zero blows straight along `+x`.
    pub yaw: f64,
    /// Pitch angle (radians) — the angle of attack, rotation of the
    /// wind toward `+z`. Positive pitches the on-coming flow upward.
    pub pitch: f64,
    /// The air the tunnel is filled with.
    pub air: Air,
    /// Upstream turbulence intensity `I` — the RMS velocity
    /// fluctuation as a fraction of the free-stream speed (`0.01` =
    /// 1 %, a smooth tunnel; `0.05`–`0.10` is on-road turbulence).
    pub turbulence_intensity: f64,
}

impl Wind {
    /// Build a straight (zero-yaw, zero-pitch) wind at the given speed
    /// in sea-level air with a 1 % turbulence intensity.
    pub fn straight(speed: f64) -> Result<Wind, AeroError> {
        Wind::new(speed, 0.0, 0.0, Air::sea_level(), 0.01)
    }

    /// Build a fully-specified wind, validating the inputs.
    pub fn new(
        speed: f64,
        yaw: f64,
        pitch: f64,
        air: Air,
        turbulence_intensity: f64,
    ) -> Result<Wind, AeroError> {
        if !(speed.is_finite() && speed >= 0.0) {
            return Err(AeroError::BadParameter {
                name: "speed",
                reason: format!("must be non-negative and finite, got {speed}"),
            });
        }
        if !yaw.is_finite() || !pitch.is_finite() {
            return Err(AeroError::IllPosedWind("yaw / pitch must be finite".into()));
        }
        if !(turbulence_intensity.is_finite() && (0.0..=1.0).contains(&turbulence_intensity)) {
            return Err(AeroError::IllPosedWind(format!(
                "turbulence intensity must be in [0, 1], got {turbulence_intensity}"
            )));
        }
        Ok(Wind {
            speed,
            yaw,
            pitch,
            air,
            turbulence_intensity,
        })
    }

    /// The unit free-stream direction vector built from the yaw /
    /// pitch angles. Zero-yaw, zero-pitch is `(1, 0, 0)`.
    pub fn direction(&self) -> Vector3<f64> {
        let (cy, sy) = (self.yaw.cos(), self.yaw.sin());
        let (cp, sp) = (self.pitch.cos(), self.pitch.sin());
        // Yaw about z, then pitch toward +z.
        Vector3::new(cy * cp, sy * cp, sp)
    }

    /// The free-stream velocity vector — `speed · direction`.
    pub fn velocity(&self) -> Vector3<f64> {
        self.speed * self.direction()
    }

    /// The dynamic pressure `q∞ = ½·ρ·U∞²` (Pa) — the normalising
    /// pressure scale for every aerodynamic coefficient.
    #[inline]
    pub fn dynamic_pressure(&self) -> f64 {
        0.5 * self.air.density * self.speed * self.speed
    }

    /// The Reynolds number for a body of characteristic length
    /// `length` — `Re = ρ·U∞·L / μ`.
    pub fn reynolds_number(&self, length: f64) -> f64 {
        if self.air.dynamic_viscosity <= 0.0 {
            return f64::INFINITY;
        }
        self.air.density * self.speed * length / self.air.dynamic_viscosity
    }

    /// The upstream turbulence kinetic energy `k∞ = 1.5·(I·U∞)²`
    /// (m²·s⁻²) — the inlet value for a turbulence model.
    pub fn inlet_tke(&self) -> f64 {
        let u_rms = self.turbulence_intensity * self.speed;
        1.5 * u_rms * u_rms
    }

    /// The upstream turbulence dissipation rate `ε∞`, estimated from
    /// `k∞` and a turbulence length scale `l` via
    /// `ε = Cμ^¾·k^{3/2}/l` with `Cμ = 0.09`.
    pub fn inlet_epsilon(&self, length_scale: f64) -> f64 {
        let k = self.inlet_tke();
        if k <= 0.0 || length_scale <= 0.0 {
            return 0.0;
        }
        let c_mu = 0.09_f64;
        c_mu.powf(0.75) * k.powf(1.5) / length_scale
    }

    /// The upstream specific dissipation rate `ω∞` for a k-ω model —
    /// `ω = k^½ / (Cμ^¼·l)`.
    pub fn inlet_omega(&self, length_scale: f64) -> f64 {
        let k = self.inlet_tke();
        if k <= 0.0 || length_scale <= 0.0 {
            return 1.0;
        }
        let c_mu = 0.09_f64;
        k.sqrt() / (c_mu.powf(0.25) * length_scale)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sea_level_air_has_standard_properties() {
        let a = Air::sea_level();
        assert!((a.density - 1.225).abs() < 1e-9);
        // ν = μ/ρ ≈ 1.48e-5 m²/s.
        let nu = a.kinematic_viscosity();
        assert!((nu - 1.81e-5 / 1.225).abs() < 1e-12);
        assert!(nu > 1.0e-5 && nu < 2.0e-5);
    }

    #[test]
    fn air_rejects_non_physical_properties() {
        assert!(Air::new(-1.0, 1e-5).is_err());
        assert!(Air::new(1.2, 0.0).is_err());
        assert!(Air::new(f64::NAN, 1e-5).is_err());
    }

    #[test]
    fn straight_wind_blows_along_x() {
        let w = Wind::straight(30.0).unwrap();
        let d = w.direction();
        assert!((d - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-12);
        assert!((w.velocity() - Vector3::new(30.0, 0.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn yaw_and_pitch_rotate_the_direction() {
        // 90° yaw points the wind along +y.
        let w = Wind::new(
            10.0,
            std::f64::consts::FRAC_PI_2,
            0.0,
            Air::sea_level(),
            0.01,
        )
        .unwrap();
        assert!((w.direction() - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-9);
        // 90° pitch points the wind straight up (+z).
        let wp = Wind::new(
            10.0,
            0.0,
            std::f64::consts::FRAC_PI_2,
            Air::sea_level(),
            0.01,
        )
        .unwrap();
        assert!((wp.direction() - Vector3::new(0.0, 0.0, 1.0)).norm() < 1e-9);
        // The direction vector is always a unit vector.
        let wd = Wind::new(20.0, 0.7, -0.3, Air::sea_level(), 0.02).unwrap();
        assert!((wd.direction().norm() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn dynamic_pressure_and_reynolds_number() {
        let w = Wind::straight(40.0).unwrap();
        // q = ½·1.225·40² = 980 Pa.
        assert!((w.dynamic_pressure() - 980.0).abs() < 1e-6);
        // Re for a 4 m car ≈ ρUL/μ — a large number, order 1e7.
        let re = w.reynolds_number(4.0);
        assert!(re > 1.0e6 && re < 2.0e7, "car Re {re} out of range");
    }

    #[test]
    fn wind_rejects_bad_speed_and_turbulence() {
        assert!(Wind::straight(-5.0).is_err());
        assert!(Wind::new(10.0, 0.0, 0.0, Air::sea_level(), 1.5).is_err());
        assert!(Wind::new(10.0, 0.0, 0.0, Air::sea_level(), -0.1).is_err());
    }

    #[test]
    fn inlet_turbulence_scalars_are_positive() {
        let w = Wind::new(50.0, 0.0, 0.0, Air::sea_level(), 0.05).unwrap();
        let k = w.inlet_tke();
        // k = 1.5·(0.05·50)² = 1.5·6.25 = 9.375.
        assert!((k - 9.375).abs() < 1e-6);
        assert!(w.inlet_epsilon(0.5) > 0.0);
        assert!(w.inlet_omega(0.5) > 0.0);
    }

    #[test]
    fn zero_speed_has_zero_dynamic_pressure() {
        let w = Wind::straight(0.0).unwrap();
        assert_eq!(w.dynamic_pressure(), 0.0);
        assert_eq!(w.inlet_tke(), 0.0);
    }
}
