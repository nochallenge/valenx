//! Defeaturing operations.

use serde::{Deserialize, Serialize};

/// One defeaturing rule.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Defeature {
    /// Detect small fillets via curvature analysis and replace each
    /// with a sharp edge.
    FilletRemove {
        /// Fillets with mean curvature radius below this are removed.
        max_radius_mm: f64,
    },
    /// Detect short cylindrical pockets and fill them.
    HoleRemove {
        /// Holes with diameter below this are filled.
        max_diameter_mm: f64,
    },
    /// Detect engraved text patches by aspect ratio and remove.
    TextRemove {
        /// Patches with depth below this are removed.
        max_depth_mm: f64,
    },
    /// Remove sliver faces (extreme aspect ratio triangles).
    SliverRemove {
        /// Triangles with `min_edge / max_edge` below this are
        /// classified as slivers.
        max_aspect: f64,
    },
}

impl Defeature {
    /// Short UI label.
    pub fn label(self) -> &'static str {
        match self {
            Self::FilletRemove { .. } => "Fillet remove",
            Self::HoleRemove { .. } => "Hole remove",
            Self::TextRemove { .. } => "Text remove",
            Self::SliverRemove { .. } => "Sliver remove",
        }
    }
}
