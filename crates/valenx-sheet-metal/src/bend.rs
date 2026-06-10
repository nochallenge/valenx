//! Bend + Flange definitions.

use serde::{Deserialize, Serialize};

/// A bend: line on the sheet (start + end in 2D outline coords) +
/// angle (rad) + inside radius (sheet units).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Bend {
    /// Bend line start `[u, v]`.
    pub start: [f64; 2],
    /// Bend line end `[u, v]`.
    pub end: [f64; 2],
    /// Bend angle in radians (positive = fold the +ve side up).
    pub angle_rad: f64,
    /// Inside radius (same units as sheet outline).
    pub inside_radius: f64,
}

impl Bend {
    /// Construct a bend.
    pub fn new(start: [f64; 2], end: [f64; 2], angle_rad: f64, inside_radius: f64) -> Self {
        Self {
            start,
            end,
            angle_rad,
            inside_radius,
        }
    }

    /// Arc length along the neutral axis for the given sheet
    /// `thickness` + `k_factor` (neutral-axis fraction, 0.0..=1.0).
    /// This is the **bend allowance** — the developed length of the
    /// arc that the unfold needs to reserve so the flat pattern
    /// folds to the correct outside dimension.
    pub fn bend_allowance(&self, thickness: f64, k_factor: f64) -> f64 {
        let r_neutral = self.inside_radius + k_factor * thickness;
        r_neutral * self.angle_rad.abs()
    }

    /// Bend deduction — the flat-pattern length removed at the bend, `BD = 2·OSSB − BA`, where the
    /// outside setback `OSSB = (inside_radius + thickness)·|tan(angle/2)|` and `BA` is the
    /// [`bend_allowance`](Self::bend_allowance). Used to size flat blanks; may be negative for
    /// shallow bends.
    pub fn bend_deduction(&self, thickness: f64, k_factor: f64) -> f64 {
        let ossb = (self.inside_radius + thickness) * (self.angle_rad / 2.0).tan().abs();
        2.0 * ossb - self.bend_allowance(thickness, k_factor)
    }
}

/// A flange: extend one outline edge by `length` at `angle_rad`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Flange {
    /// Which outline edge to attach to (0 = first segment between
    /// `outline[0]` and `outline[1]`, wrapping around at the last
    /// vertex).
    pub edge_id: usize,
    /// Flange length (sheet units).
    pub length: f64,
    /// Bend angle at the root of the flange (rad).
    pub angle_rad: f64,
}

impl Flange {
    /// Construct a flange.
    pub fn new(edge_id: usize, length: f64, angle_rad: f64) -> Self {
        Self {
            edge_id,
            length,
            angle_rad,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bend_allowance_90deg() {
        let b = Bend::new([0.0, 0.0], [1.0, 0.0], std::f64::consts::FRAC_PI_2, 1.0);
        let ba = b.bend_allowance(1.0, 0.44);
        // 90° arc of neutral-axis radius 1.44 → 1.44 * π/2.
        assert!((ba - 1.44 * std::f64::consts::FRAC_PI_2).abs() < 1e-9);
    }

    #[test]
    fn bend_deduction_90deg() {
        // 90° bend, r=5, t=2, K=0.44: OSSB=(5+2)·tan(45°)=7, BA=5.88·π/2≈9.2363, BD≈4.7637.
        let b = Bend::new([0.0, 0.0], [1.0, 0.0], std::f64::consts::FRAC_PI_2, 5.0);
        let bd = b.bend_deduction(2.0, 0.44);
        assert!((bd - 4.7637).abs() < 1e-3);
        // BD = 2·OSSB − BA (composes the shipped bend_allowance).
        let ossb = 7.0 * std::f64::consts::FRAC_PI_4.tan();
        assert!((bd - (2.0 * ossb - b.bend_allowance(2.0, 0.44))).abs() < 1e-9);
    }
}
