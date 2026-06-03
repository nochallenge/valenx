//! Piping system — network of sections + junction-mounted
//! fittings / valves.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::fitting::{PipeFitting, Valve};
use crate::pipe::PipeSection;

/// A fitting or valve placed at a specific world-space point in the
/// network.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Junction {
    /// World-space location.
    pub at: Vector3<f64>,
    /// One of `Fitting`, `ValveAt`.
    pub kind: JunctionKind,
}

/// Discriminant for [`Junction::kind`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum JunctionKind {
    /// Plain pipe fitting.
    Fitting(PipeFitting),
    /// Inline valve.
    ValveAt(Valve),
}

/// A piping system — every section + every junction.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Piping {
    /// All pipe sections in the network.
    pub sections: Vec<PipeSection>,
    /// All fittings / valves placed at junctions.
    pub junctions: Vec<Junction>,
}

impl Piping {
    /// Empty network.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a section.
    pub fn push_section(&mut self, s: PipeSection) {
        self.sections.push(s);
    }

    /// Append a junction.
    pub fn push_junction(&mut self, j: Junction) {
        self.junctions.push(j);
    }

    /// Total piping run length in mm.
    pub fn total_length_mm(&self) -> f64 {
        self.sections.iter().map(|s| s.length_mm()).sum()
    }

    /// Count of fittings and valves combined.
    pub fn fitting_count(&self) -> usize {
        self.junctions.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dims::Schedule;
    use crate::pipe::Material;

    fn demo_section() -> PipeSection {
        PipeSection::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1000.0),
            "2",
            Schedule::Sch40,
            Material::CarbonSteel,
        )
    }

    #[test]
    fn total_length_sums_sections() {
        let mut p = Piping::new();
        p.push_section(demo_section());
        p.push_section(demo_section());
        assert!((p.total_length_mm() - 2000.0).abs() < 1e-6);
    }

    #[test]
    fn fitting_count_tracks_junctions() {
        let mut p = Piping::new();
        p.push_junction(Junction {
            at: Vector3::zeros(),
            kind: JunctionKind::Fitting(PipeFitting::Elbow90),
        });
        p.push_junction(Junction {
            at: Vector3::new(100.0, 0.0, 0.0),
            kind: JunctionKind::ValveAt(Valve::Ball),
        });
        assert_eq!(p.fitting_count(), 2);
    }
}
