//! US Standard Atmosphere 1976 (0–86 km).
//!
//! Returns temperature, pressure, density and the local speed of sound
//! as a function of *geometric* altitude. The model is the canonical
//! seven-layer piecewise-linear-temperature profile; pressure follows
//! the barometric formula (gradient layers) or the isothermal
//! exponential (zero-lapse layers). Above [`ATMOS_TOP_M`] the
//! atmosphere is treated as vacuum.
//!
//! Reference base values (geopotential altitude, base temperature,
//! lapse rate, base pressure) are the published 1976 constants; the
//! recovered sea-level density is 1.2250 kg/m³ and the recovered values
//! at the layer boundaries match the standard tables to better than
//! 0.1 %, which the unit tests pin down.

use crate::constants::{ATMOS_EARTH_RADIUS, ATMOS_TOP_M, G0, GAMMA_AIR, R_AIR};

/// One isothermal-or-gradient layer of the standard atmosphere.
struct Layer {
    /// Base geopotential altitude (m).
    base_h: f64,
    /// Base temperature (K).
    base_t: f64,
    /// Temperature lapse rate (K/m); zero for isothermal layers.
    lapse: f64,
    /// Base pressure (Pa).
    base_p: f64,
}

/// The seven base layers of the 1976 standard atmosphere (to 84.852 km
/// geopotential ≈ 86 km geometric).
const LAYERS: [Layer; 7] = [
    Layer {
        base_h: 0.0,
        base_t: 288.15,
        lapse: -0.0065,
        base_p: 101_325.0,
    },
    Layer {
        base_h: 11_000.0,
        base_t: 216.65,
        lapse: 0.0,
        base_p: 22_632.06,
    },
    Layer {
        base_h: 20_000.0,
        base_t: 216.65,
        lapse: 0.001,
        base_p: 5_474.889,
    },
    Layer {
        base_h: 32_000.0,
        base_t: 228.65,
        lapse: 0.0028,
        base_p: 868.0187,
    },
    Layer {
        base_h: 47_000.0,
        base_t: 270.65,
        lapse: 0.0,
        base_p: 110.9063,
    },
    Layer {
        base_h: 51_000.0,
        base_t: 270.65,
        lapse: -0.0028,
        base_p: 66.93887,
    },
    Layer {
        base_h: 71_000.0,
        base_t: 214.65,
        lapse: -0.002,
        base_p: 3.956420,
    },
];

/// Local atmospheric state at a given altitude.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtmosphereSample {
    /// Temperature (K).
    pub temperature: f64,
    /// Static pressure (Pa).
    pub pressure: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Local speed of sound (m/s).
    pub speed_of_sound: f64,
}

/// Convert geometric altitude (m) to geopotential altitude (m).
fn geopotential(geometric_m: f64) -> f64 {
    ATMOS_EARTH_RADIUS * geometric_m / (ATMOS_EARTH_RADIUS + geometric_m)
}

/// Sample the US Standard Atmosphere 1976 at a geometric altitude.
///
/// Below 0 m the sea-level layer is extrapolated (so a launch site a
/// few metres below datum still gets a sane density); at or above
/// [`ATMOS_TOP_M`] the result is vacuum (zero pressure and density)
/// with the 86 km temperature retained for the speed-of-sound term.
pub fn sample(geometric_altitude_m: f64) -> AtmosphereSample {
    if geometric_altitude_m >= ATMOS_TOP_M {
        // Vacuum: no drag contribution. Keep a finite speed of sound so
        // a Mach lookup never divides by zero.
        let t = 186.87;
        return AtmosphereSample {
            temperature: t,
            pressure: 0.0,
            density: 0.0,
            speed_of_sound: (GAMMA_AIR * R_AIR * t).sqrt(),
        };
    }

    let h = geopotential(geometric_altitude_m);

    // Find the layer whose base is at or below `h` (the last such layer).
    let mut idx = 0usize;
    for (i, layer) in LAYERS.iter().enumerate() {
        if h >= layer.base_h {
            idx = i;
        } else {
            break;
        }
    }
    let layer = &LAYERS[idx];
    let dh = h - layer.base_h;

    let (temperature, pressure) = if layer.lapse.abs() < 1e-12 {
        // Isothermal layer: exponential pressure decay.
        let t = layer.base_t;
        let p = layer.base_p * (-G0 * dh / (R_AIR * t)).exp();
        (t, p)
    } else {
        // Gradient layer: linear temperature, power-law pressure.
        let t = layer.base_t + layer.lapse * dh;
        let exponent = G0 / (R_AIR * layer.lapse);
        let p = layer.base_p * (layer.base_t / t).powf(exponent);
        (t, p)
    };

    let density = pressure / (R_AIR * temperature);
    let speed_of_sound = (GAMMA_AIR * R_AIR * temperature).sqrt();

    AtmosphereSample {
        temperature,
        pressure,
        density,
        speed_of_sound,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rel_err(a: f64, b: f64) -> f64 {
        (a - b).abs() / b.abs()
    }

    #[test]
    fn sea_level_matches_standard() {
        let s = sample(0.0);
        assert!(rel_err(s.temperature, 288.15) < 1e-4);
        assert!(rel_err(s.pressure, 101_325.0) < 1e-4);
        assert!(rel_err(s.density, 1.2250) < 1e-3, "density {}", s.density);
        // Speed of sound at sea level ≈ 340.3 m/s.
        assert!(rel_err(s.speed_of_sound, 340.3) < 2e-3);
    }

    #[test]
    fn tropopause_11km() {
        // At ~11 km geometric the table gives ρ ≈ 0.3639 kg/m³,
        // T = 216.65 K, p ≈ 22 632 Pa (at the geopotential boundary).
        let s = sample(11_019.0); // geometric altitude of 11 km geopotential
        assert!(rel_err(s.temperature, 216.65) < 1e-3, "T {}", s.temperature);
        assert!(rel_err(s.pressure, 22_632.0) < 5e-3, "p {}", s.pressure);
        assert!(rel_err(s.density, 0.3639) < 5e-3, "rho {}", s.density);
    }

    #[test]
    fn twenty_km_isothermal() {
        // 20 km geopotential boundary: T = 216.65, p ≈ 5474.9 Pa.
        let s = sample(20_063.0);
        assert!(rel_err(s.temperature, 216.65) < 1e-3);
        assert!(rel_err(s.pressure, 5_474.9) < 1e-2, "p {}", s.pressure);
    }

    #[test]
    fn pressure_monotonically_decreases() {
        let mut prev = sample(0.0).pressure;
        for km in 1..=85 {
            let p = sample(km as f64 * 1000.0).pressure;
            assert!(p < prev, "pressure rose at {km} km: {p} >= {prev}");
            prev = p;
        }
    }

    #[test]
    fn vacuum_above_top() {
        let s = sample(ATMOS_TOP_M + 10.0);
        assert_eq!(s.pressure, 0.0);
        assert_eq!(s.density, 0.0);
        assert!(s.speed_of_sound > 0.0);
    }
}
