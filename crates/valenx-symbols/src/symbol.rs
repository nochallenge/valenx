//! Symbol kind library — canonical electrical, hydraulic, and
//! pneumatic glyphs as SVG path strings. Each glyph is centred on
//! the origin and sized roughly 1 logical unit across (the caller
//! scales / translates via [`crate::schematic::PlacedSymbol`]).

use serde::{Deserialize, Serialize};

/// Coarse symbol family.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SymbolFamily {
    /// Electrical schematic.
    Electrical,
    /// Hydraulic / fluid power.
    Hydraulic,
    /// Pneumatic / gas power.
    Pneumatic,
}

/// All supported schematic glyphs.
///
/// Three families × 20+ total entries covers the standard FreeCAD
/// `Symbols Library` community workbench surface. Each variant has a
/// canonical SVG path string emitted by [`SymbolKind::to_svg_path`]
/// and a human-readable [`SymbolKind::label`] for the UI.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SymbolKind {
    // ---- Electrical ----
    /// Zig-zag resistor.
    Resistor,
    /// Parallel-plate capacitor.
    Capacitor,
    /// Coil inductor.
    Inductor,
    /// Triangle + bar diode.
    Diode,
    /// Light-emitting diode (Diode + arrows).
    Led,
    /// Bipolar transistor.
    Transistor,
    /// Single-throw switch.
    Switch,
    /// Ground reference.
    Ground,
    /// Positive supply.
    VPlus,
    /// Battery cell pair.
    Battery,
    /// Circular motor.
    Motor,
    /// Incandescent lamp.
    Lamp,

    // ---- Hydraulic ----
    /// Hydraulic pump (circle + triangle arrow).
    HydraulicPump,
    /// Spool valve (rectangle with arrows).
    HydraulicValve,
    /// Single-acting cylinder.
    HydraulicCylinder,
    /// Open-top reservoir.
    HydraulicReservoir,
    /// Diamond-shape filter.
    HydraulicFilter,

    // ---- Pneumatic ----
    /// Pneumatic compressor (circle + triangle pointing out).
    PneumaticCompressor,
    /// Pneumatic cylinder (double-acting).
    PneumaticCylinder,
    /// Pneumatic 5/2 valve.
    PneumaticValve,
    /// Pressure regulator (square + arrow).
    PneumaticRegulator,
}

impl SymbolKind {
    /// Symbol family.
    pub fn family(self) -> SymbolFamily {
        match self {
            Self::Resistor
            | Self::Capacitor
            | Self::Inductor
            | Self::Diode
            | Self::Led
            | Self::Transistor
            | Self::Switch
            | Self::Ground
            | Self::VPlus
            | Self::Battery
            | Self::Motor
            | Self::Lamp => SymbolFamily::Electrical,
            Self::HydraulicPump
            | Self::HydraulicValve
            | Self::HydraulicCylinder
            | Self::HydraulicReservoir
            | Self::HydraulicFilter => SymbolFamily::Hydraulic,
            Self::PneumaticCompressor
            | Self::PneumaticCylinder
            | Self::PneumaticValve
            | Self::PneumaticRegulator => SymbolFamily::Pneumatic,
        }
    }

    /// Short UI label.
    pub fn label(self) -> &'static str {
        match self {
            Self::Resistor => "Resistor",
            Self::Capacitor => "Capacitor",
            Self::Inductor => "Inductor",
            Self::Diode => "Diode",
            Self::Led => "LED",
            Self::Transistor => "Transistor",
            Self::Switch => "Switch",
            Self::Ground => "Ground (GND)",
            Self::VPlus => "V+",
            Self::Battery => "Battery",
            Self::Motor => "Motor",
            Self::Lamp => "Lamp",
            Self::HydraulicPump => "Hydraulic pump",
            Self::HydraulicValve => "Hydraulic valve",
            Self::HydraulicCylinder => "Hydraulic cylinder",
            Self::HydraulicReservoir => "Hydraulic reservoir",
            Self::HydraulicFilter => "Hydraulic filter",
            Self::PneumaticCompressor => "Pneumatic compressor",
            Self::PneumaticCylinder => "Pneumatic cylinder",
            Self::PneumaticValve => "Pneumatic 5/2 valve",
            Self::PneumaticRegulator => "Pneumatic regulator",
        }
    }

    /// Canonical SVG path string. Paths are designed for a 60-unit
    /// bounding box centred on (0,0). The schematic renderer wraps
    /// each glyph in a `<g transform="translate(...) rotate(...)">`.
    pub fn to_svg_path(self) -> &'static str {
        match self {
            // Electrical
            Self::Resistor => {
                // Zig-zag with two pigtail leads.
                "M -30 0 L -20 0 L -16 -8 L -8 8 L 0 -8 L 8 8 L 16 -8 L 20 0 L 30 0"
            }
            Self::Capacitor => {
                // Two parallel plates + leads.
                "M -30 0 L -4 0 M -4 -10 L -4 10 M 4 -10 L 4 10 M 4 0 L 30 0"
            }
            Self::Inductor => {
                // Four loops + leads.
                "M -30 0 L -20 0 \
                 a 5 5 0 0 1 10 0 \
                 a 5 5 0 0 1 10 0 \
                 a 5 5 0 0 1 10 0 \
                 a 5 5 0 0 1 10 0 \
                 L 30 0"
            }
            Self::Diode => {
                // Triangle pointing right + cathode bar.
                "M -30 0 L -10 0 L -10 -10 L 10 0 L -10 10 L -10 0 M 10 -10 L 10 10 M 10 0 L 30 0"
            }
            Self::Led => {
                // Diode + two arrows.
                "M -30 0 L -10 0 L -10 -10 L 10 0 L -10 10 L -10 0 M 10 -10 L 10 10 M 10 0 L 30 0 \
                 M 16 -16 L 24 -24 M 22 -24 L 24 -24 L 24 -22 \
                 M 8 -16 L 16 -24 M 14 -24 L 16 -24 L 16 -22"
            }
            Self::Transistor => {
                // NPN bipolar: circle + base/collector/emitter.
                "M 0 0 m -15 0 a 15 15 0 1 0 30 0 a 15 15 0 1 0 -30 0 \
                 M -30 0 L -8 0 \
                 M -8 -10 L -8 10 \
                 M -8 -3 L 12 -15 L 12 -30 \
                 M -8 3 L 12 15 L 12 30 \
                 M 6 9 L 12 15 L 6 15"
            }
            Self::Switch => {
                // Open SPST: dot + tilted bar + dot.
                "M -30 0 L -10 0 M 10 0 L 30 0 \
                 M -10 0 L 12 -12 \
                 M -10 0 m -2 0 a 2 2 0 1 0 4 0 a 2 2 0 1 0 -4 0 \
                 M 10 0 m -2 0 a 2 2 0 1 0 4 0 a 2 2 0 1 0 -4 0"
            }
            Self::Ground => {
                // Three horizontal lines getting shorter.
                "M 0 -20 L 0 0 \
                 M -15 0 L 15 0 \
                 M -10 6 L 10 6 \
                 M -5 12 L 5 12"
            }
            Self::VPlus => {
                // Vertical line with a + sign above.
                "M 0 0 L 0 -20 \
                 M -6 -25 L 6 -25 \
                 M 0 -31 L 0 -19"
            }
            Self::Battery => {
                // Two-cell pair: long-short-long-short.
                "M -30 0 L -10 0 \
                 M -10 -12 L -10 12 \
                 M -5 -6 L -5 6 \
                 M 5 -12 L 5 12 \
                 M 10 -6 L 10 6 \
                 M 10 0 L 30 0"
            }
            Self::Motor => {
                // Circle with M inside.
                "M 0 0 m -15 0 a 15 15 0 1 0 30 0 a 15 15 0 1 0 -30 0 \
                 M -30 0 L -15 0 \
                 M 15 0 L 30 0 \
                 M -8 7 L -8 -7 L 0 0 L 8 -7 L 8 7"
            }
            Self::Lamp => {
                // Circle with X inside.
                "M 0 0 m -15 0 a 15 15 0 1 0 30 0 a 15 15 0 1 0 -30 0 \
                 M -30 0 L -15 0 \
                 M 15 0 L 30 0 \
                 M -10 -10 L 10 10 \
                 M -10 10 L 10 -10"
            }

            // Hydraulic
            Self::HydraulicPump => {
                // Circle + filled triangle pointing right.
                "M 0 0 m -20 0 a 20 20 0 1 0 40 0 a 20 20 0 1 0 -40 0 \
                 M -10 -8 L 10 0 L -10 8 Z"
            }
            Self::HydraulicValve => {
                // Two-rectangle 4/3 spool valve.
                "M -30 -15 L 0 -15 L 0 15 L -30 15 Z \
                 M 0 -15 L 30 -15 L 30 15 L 0 15 \
                 M -25 0 L -5 0 \
                 M 5 0 L 25 0 \
                 M -25 -10 L -10 -10 L -10 10 L -25 10"
            }
            Self::HydraulicCylinder => {
                // Long rectangle with piston + rod.
                "M -30 -10 L 30 -10 L 30 10 L -30 10 Z \
                 M 10 -10 L 10 10 \
                 M 10 0 L 35 0"
            }
            Self::HydraulicReservoir => {
                // Open-top trapezoid.
                "M -25 -15 L -25 15 L 25 15 L 25 -15 \
                 M -25 0 L -35 0"
            }
            Self::HydraulicFilter => {
                // Diamond with dashed centre.
                "M 0 -20 L 20 0 L 0 20 L -20 0 Z \
                 M -10 0 L 10 0"
            }

            // Pneumatic
            Self::PneumaticCompressor => {
                // Circle + triangle pointing outward (left).
                "M 0 0 m -20 0 a 20 20 0 1 0 40 0 a 20 20 0 1 0 -40 0 \
                 M 10 -8 L -10 0 L 10 8 Z"
            }
            Self::PneumaticCylinder => {
                // Double-acting cylinder.
                "M -30 -10 L 30 -10 L 30 10 L -30 10 Z \
                 M 0 -10 L 0 10 \
                 M -30 -3 L -35 -3 M -30 3 L -35 3 \
                 M 0 0 L 35 0"
            }
            Self::PneumaticValve => {
                // 5/2 directional valve (2 boxes side by side).
                "M -30 -15 L 0 -15 L 0 15 L -30 15 Z \
                 M 0 -15 L 30 -15 L 30 15 L 0 15 \
                 M -25 -5 L -5 -5 \
                 M -25 5 L -5 5 \
                 M 5 -5 L 25 5 \
                 M 5 5 L 25 -5"
            }
            Self::PneumaticRegulator => {
                // Square + slanted arrow.
                "M -15 -15 L 15 -15 L 15 15 L -15 15 Z \
                 M -10 10 L 10 -10 \
                 M 6 -10 L 10 -10 L 10 -6"
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn library_covers_20_plus_kinds() {
        let kinds = [
            SymbolKind::Resistor,
            SymbolKind::Capacitor,
            SymbolKind::Inductor,
            SymbolKind::Diode,
            SymbolKind::Led,
            SymbolKind::Transistor,
            SymbolKind::Switch,
            SymbolKind::Ground,
            SymbolKind::VPlus,
            SymbolKind::Battery,
            SymbolKind::Motor,
            SymbolKind::Lamp,
            SymbolKind::HydraulicPump,
            SymbolKind::HydraulicValve,
            SymbolKind::HydraulicCylinder,
            SymbolKind::HydraulicReservoir,
            SymbolKind::HydraulicFilter,
            SymbolKind::PneumaticCompressor,
            SymbolKind::PneumaticCylinder,
            SymbolKind::PneumaticValve,
            SymbolKind::PneumaticRegulator,
        ];
        assert_eq!(kinds.len(), 21);
        for k in kinds {
            assert!(!k.to_svg_path().is_empty(), "{k:?} has no SVG path");
            assert!(!k.label().is_empty(), "{k:?} has no label");
        }
    }

    #[test]
    fn families_partition_kinds() {
        assert_eq!(SymbolKind::Resistor.family(), SymbolFamily::Electrical);
        assert_eq!(SymbolKind::HydraulicPump.family(), SymbolFamily::Hydraulic);
        assert_eq!(SymbolKind::PneumaticValve.family(), SymbolFamily::Pneumatic);
    }
}
