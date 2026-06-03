//! `KicadBoard` — parametric description of a KiCad PCB.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Pad geometry.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum PadShape {
    /// Round pad.
    Circle,
    /// Rectangular pad.
    Rect,
    /// Oval / stadium pad (size_x != size_y, rounded ends).
    Oval,
}

/// One pad on the PCB.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Pad {
    /// Centre position in board coordinates (mm).
    pub position: Vector3<f64>,
    /// Pad shape.
    pub shape: PadShape,
    /// (size_x, size_y) in mm.
    pub size_mm: [f64; 2],
    /// Layer name (e.g. "F.Cu", "B.Cu").
    pub layer: String,
}

/// One placed component (footprint instance).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Component {
    /// Reference designator (e.g. "R1", "U2").
    pub ref_designator: String,
    /// Footprint name (e.g. "Resistor_SMD:R_0805").
    pub footprint_name: String,
    /// Centre position in board coordinates (mm).
    pub position: Vector3<f64>,
    /// Z-axis rotation in degrees.
    pub rotation_deg: f64,
    /// Path to the .step / .wrl 3D model file (optional).
    pub model_3d_path: Option<String>,
}

/// One PCB.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct KicadBoard {
    /// 2D outline polygon (closed CCW, mm).
    pub outline: Vec<[f64; 2]>,
    /// Board thickness (mm). Default 1.6 mm.
    pub thickness_mm: f64,
    /// Through-hole drills: (position, diameter_mm).
    pub drill_holes: Vec<(Vector3<f64>, f64)>,
    /// Surface-mount + through-hole pads.
    pub pads: Vec<Pad>,
    /// Placed components.
    pub components: Vec<Component>,
}

impl KicadBoard {
    /// Empty board with the standard 1.6 mm FR-4 thickness.
    pub fn new_default() -> Self {
        Self {
            outline: Vec::new(),
            thickness_mm: 1.6,
            drill_holes: Vec::new(),
            pads: Vec::new(),
            components: Vec::new(),
        }
    }

    /// A simple 100×80 mm rectangular dev board with a corner drill.
    pub fn demo_devboard() -> Self {
        let mut b = Self::new_default();
        b.outline = vec![
            [0.0, 0.0],
            [100.0, 0.0],
            [100.0, 80.0],
            [0.0, 80.0],
        ];
        b.drill_holes.push((Vector3::new(5.0, 5.0, 0.0), 3.2));
        b.drill_holes.push((Vector3::new(95.0, 5.0, 0.0), 3.2));
        b.drill_holes.push((Vector3::new(5.0, 75.0, 0.0), 3.2));
        b.drill_holes.push((Vector3::new(95.0, 75.0, 0.0), 3.2));
        b.components.push(Component {
            ref_designator: "U1".into(),
            footprint_name: "Package_DIP:DIP-8_W7.62mm".into(),
            position: Vector3::new(50.0, 40.0, 0.0),
            rotation_deg: 0.0,
            model_3d_path: Some("DIP-8.step".into()),
        });
        b
    }
}
