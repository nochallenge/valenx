//! HVAC equipment kinds + parametric CAD solids.

use serde::{Deserialize, Serialize};

use valenx_cad::primitives::{box_solid, cylinder};
use valenx_cad::Solid;

use crate::error::HvacError;

/// HVAC equipment kind.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Equipment {
    /// Air Handling Unit.
    Ahu,
    /// Variable Air Volume box.
    Vav,
    /// Supply diffuser.
    Diffuser,
    /// Return grille.
    Grille,
    /// Fan.
    Fan,
    /// Heater (electric / hydronic).
    Heater,
    /// Chiller.
    Chiller,
    /// Damper.
    Damper,
}

impl Equipment {
    /// Short label for UI display.
    pub fn label(self) -> &'static str {
        match self {
            Equipment::Ahu => "AHU",
            Equipment::Vav => "VAV",
            Equipment::Diffuser => "Diffuser",
            Equipment::Grille => "Grille",
            Equipment::Fan => "Fan",
            Equipment::Heater => "Heater",
            Equipment::Chiller => "Chiller",
            Equipment::Damper => "Damper",
        }
    }

    /// Default bounding-box size in mm (w, h, l) for the equipment
    /// kind. Used when callers don't provide one explicitly.
    pub fn default_size_mm(self) -> (f64, f64, f64) {
        match self {
            Equipment::Ahu => (1500.0, 1500.0, 2000.0),
            Equipment::Vav => (400.0, 400.0, 600.0),
            Equipment::Diffuser => (300.0, 100.0, 300.0),
            Equipment::Grille => (300.0, 300.0, 50.0),
            Equipment::Fan => (500.0, 500.0, 500.0),
            Equipment::Heater => (400.0, 600.0, 800.0),
            Equipment::Chiller => (2000.0, 1500.0, 3000.0),
            Equipment::Damper => (300.0, 300.0, 100.0),
        }
    }
}

/// Emit a CAD [`Solid`] for the given equipment kind. v1: rectangular
/// box for prismatic equipment, cylinder for the round Fan / Damper
/// types. `(w, h, l)` are in mm.
pub fn to_solid(kind: Equipment, size: (f64, f64, f64)) -> Result<Solid, HvacError> {
    let (w, h, l) = size;
    if w <= 0.0 || h <= 0.0 || l <= 0.0 {
        return Err(HvacError::BadParameter {
            name: "size",
            reason: format!("all dimensions must be > 0, got ({w}, {h}, {l})"),
        });
    }
    match kind {
        Equipment::Fan | Equipment::Damper => {
            cylinder(w / 2.0, l).map_err(|e| HvacError::Cad(e.to_string()))
        }
        _ => box_solid(w, h, l).map_err(|e| HvacError::Cad(e.to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn label_is_non_empty_for_every_kind() {
        for k in [
            Equipment::Ahu,
            Equipment::Vav,
            Equipment::Diffuser,
            Equipment::Grille,
            Equipment::Fan,
            Equipment::Heater,
            Equipment::Chiller,
            Equipment::Damper,
        ] {
            assert!(!k.label().is_empty());
        }
    }

    #[test]
    fn default_sizes_are_positive() {
        for k in [Equipment::Ahu, Equipment::Vav, Equipment::Fan] {
            let (w, h, l) = k.default_size_mm();
            assert!(w > 0.0 && h > 0.0 && l > 0.0);
        }
    }

    #[test]
    fn to_solid_emits_cylinder_for_fan() {
        let s = to_solid(Equipment::Fan, (300.0, 300.0, 200.0)).unwrap();
        assert!(s.faces() > 0);
    }

    #[test]
    fn to_solid_rejects_zero_dimension() {
        let err = to_solid(Equipment::Ahu, (0.0, 100.0, 100.0)).unwrap_err();
        assert!(matches!(err, HvacError::BadParameter { .. }));
    }
}
