//! Light sources.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// One light. All colours are linear RGB in `[0, 1]`; intensity is in
/// watts.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Light {
    /// Omnidirectional point light.
    Point {
        /// World position.
        position: Vector3<f64>,
        /// Linear RGB.
        color: [f32; 3],
        /// Watts.
        intensity: f32,
    },
    /// Parallel (sun-like) light from a direction.
    Directional {
        /// Light direction (points *towards* the scene).
        direction: Vector3<f64>,
        /// Linear RGB.
        color: [f32; 3],
        /// Watts / m².
        irradiance: f32,
    },
    /// Cone-shaped spot light.
    Spot {
        /// World position.
        position: Vector3<f64>,
        /// Light direction (axis of the cone).
        direction: Vector3<f64>,
        /// Half-angle of the inner cone in radians.
        inner_angle_rad: f64,
        /// Half-angle of the outer cone (falloff) in radians.
        outer_angle_rad: f64,
        /// Linear RGB.
        color: [f32; 3],
        /// Watts.
        intensity: f32,
    },
    /// Rectangular area emitter.
    Area {
        /// Corner.
        corner: Vector3<f64>,
        /// First edge from `corner`.
        edge_u: Vector3<f64>,
        /// Second edge from `corner`.
        edge_v: Vector3<f64>,
        /// Linear RGB.
        color: [f32; 3],
        /// Watts / m².
        radiance: f32,
    },
}

impl Light {
    /// Short label for UI dropdowns + error messages.
    pub fn label(&self) -> &'static str {
        match self {
            Self::Point { .. } => "Point",
            Self::Directional { .. } => "Directional",
            Self::Spot { .. } => "Spot",
            Self::Area { .. } => "Area",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one_of_each() -> Vec<Light> {
        vec![
            Light::Point {
                position: Vector3::new(0.0, 5.0, 0.0),
                color: [1.0, 1.0, 1.0],
                intensity: 100.0,
            },
            Light::Directional {
                direction: Vector3::new(0.0, -1.0, 0.0),
                color: [1.0, 0.95, 0.9],
                irradiance: 5.0,
            },
            Light::Spot {
                position: Vector3::new(2.0, 4.0, 1.0),
                direction: Vector3::new(0.0, -1.0, 0.0),
                inner_angle_rad: 0.3,
                outer_angle_rad: 0.5,
                color: [1.0, 1.0, 0.8],
                intensity: 60.0,
            },
            Light::Area {
                corner: Vector3::new(-1.0, 3.0, -1.0),
                edge_u: Vector3::new(2.0, 0.0, 0.0),
                edge_v: Vector3::new(0.0, 0.0, 2.0),
                color: [1.0, 1.0, 1.0],
                radiance: 10.0,
            },
        ]
    }

    #[test]
    fn every_variant_has_a_distinct_label() {
        // Drives all four arms of `Light::label`.
        let labels: Vec<&str> = one_of_each().iter().map(|l| l.label()).collect();
        assert_eq!(labels, ["Point", "Directional", "Spot", "Area"]);
    }

    #[test]
    fn lights_round_trip_through_ron() {
        for light in one_of_each() {
            let ron = ron::to_string(&light).unwrap();
            let back: Light = ron::from_str(&ron).unwrap();
            assert_eq!(light, back, "{} should round-trip", light.label());
        }
    }

    #[test]
    fn light_is_clone_and_comparable() {
        let lights = one_of_each();
        let point = lights[0].clone();
        assert_eq!(point, lights[0]);
        // Different variants are never equal.
        assert_ne!(lights[0], lights[1]);
        // Debug formatting works (used in failure messages).
        assert!(format!("{point:?}").contains("Point"));
    }
}
