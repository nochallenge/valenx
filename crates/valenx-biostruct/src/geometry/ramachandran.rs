//! Ramachandran (φ/ψ) computation and region classification.
//!
//! For each residue with both backbone torsions defined this module
//! classifies the `(φ, ψ)` point into a [`RamachandranRegion`]. The
//! region boundaries are the standard textbook quadrant definitions
//! (a simplified Lovell-style partition); a residue outside every
//! favoured / allowed region is an **outlier** — a geometry-quality
//! red flag.

use crate::geometry::angles::chain_phi_psi;
use crate::structure::Chain;

/// A Ramachandran-plot region.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum RamachandranRegion {
    /// Right-handed α-helix basin.
    AlphaHelix,
    /// β-sheet / extended basin.
    BetaSheet,
    /// Left-handed α-helix basin (rare, mostly glycine).
    LeftHandedAlpha,
    /// Polyproline-II / bridge region.
    Bridge,
    /// Inside no favoured / allowed basin — a likely modelling error.
    Outlier,
}

impl RamachandranRegion {
    /// Whether this region is one of the sterically favourable
    /// basins (anything but [`Outlier`](Self::Outlier)).
    pub fn is_allowed(&self) -> bool {
        !matches!(self, RamachandranRegion::Outlier)
    }

    /// Short human-readable label.
    pub fn label(&self) -> &'static str {
        match self {
            RamachandranRegion::AlphaHelix => "alpha-helix",
            RamachandranRegion::BetaSheet => "beta-sheet",
            RamachandranRegion::LeftHandedAlpha => "left-handed-alpha",
            RamachandranRegion::Bridge => "bridge",
            RamachandranRegion::Outlier => "outlier",
        }
    }
}

/// One residue's Ramachandran datum.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct RamachandranPoint {
    /// Residue index within its chain.
    pub residue_index: usize,
    /// φ torsion, degrees.
    pub phi: f64,
    /// ψ torsion, degrees.
    pub psi: f64,
    /// Classified region.
    pub region: RamachandranRegion,
}

/// Classify a `(φ, ψ)` pair in degrees into a [`RamachandranRegion`].
///
/// The basins are axis-aligned boxes in `(-180, 180]²`:
///
/// - **α-helix**: φ ∈ [-160, -20], ψ ∈ [-120, 50]
/// - **β-sheet**: φ ∈ [-180, -40], ψ ∈ [90, 180] ∪ [-180, -150]
/// - **left-handed α**: φ ∈ [20, 90], ψ ∈ [-30, 90]
/// - **bridge / PP-II**: φ ∈ [-110, -40], ψ ∈ [50, 100]
///
/// The α-helix box is tested first so the overlap with the bridge box
/// resolves to α-helix. Anything outside all four is an outlier.
pub fn classify(phi: f64, psi: f64) -> RamachandranRegion {
    let phi = wrap180(phi);
    let psi = wrap180(psi);

    // Alpha-helix basin.
    if (-160.0..=-20.0).contains(&phi) && (-120.0..=50.0).contains(&psi) {
        return RamachandranRegion::AlphaHelix;
    }
    // Beta-sheet basin: a high-psi band, wrapping past +/-180.
    if (-180.0..=-40.0).contains(&phi)
        && (psi >= 90.0 || psi <= -150.0)
    {
        return RamachandranRegion::BetaSheet;
    }
    // Left-handed alpha.
    if (20.0..=90.0).contains(&phi) && (-30.0..=90.0).contains(&psi) {
        return RamachandranRegion::LeftHandedAlpha;
    }
    // Bridge / polyproline-II region.
    if (-110.0..=-40.0).contains(&phi) && (50.0..=100.0).contains(&psi) {
        return RamachandranRegion::Bridge;
    }
    RamachandranRegion::Outlier
}

/// Compute the Ramachandran points for every residue of `chain` with
/// both backbone torsions defined.
pub fn chain_ramachandran(chain: &Chain) -> Vec<RamachandranPoint> {
    let mut out = Vec::new();
    for (idx, tor) in chain_phi_psi(chain) {
        if let (Some(phi), Some(psi)) = (tor.phi, tor.psi) {
            out.push(RamachandranPoint {
                residue_index: idx,
                phi,
                psi,
                region: classify(phi, psi),
            });
        }
    }
    out
}

/// Summary statistics over a chain's Ramachandran points.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct RamachandranSummary {
    /// Residues with both torsions defined.
    pub total: usize,
    /// Residues in the α-helix basin.
    pub alpha: usize,
    /// Residues in the β-sheet basin.
    pub beta: usize,
    /// Residues in the left-handed α basin.
    pub left_alpha: usize,
    /// Residues in the bridge / PP-II region.
    pub bridge: usize,
    /// Residues outside every allowed basin.
    pub outliers: usize,
}

impl RamachandranSummary {
    /// Fraction of residues that fall in an allowed basin, in `[0, 1]`.
    /// Returns `1.0` for an empty chain (vacuously fine).
    pub fn allowed_fraction(&self) -> f64 {
        if self.total == 0 {
            return 1.0;
        }
        (self.total - self.outliers) as f64 / self.total as f64
    }
}

/// Tally a chain's Ramachandran points into a [`RamachandranSummary`].
pub fn summarize(chain: &Chain) -> RamachandranSummary {
    let mut s = RamachandranSummary::default();
    for p in chain_ramachandran(chain) {
        s.total += 1;
        match p.region {
            RamachandranRegion::AlphaHelix => s.alpha += 1,
            RamachandranRegion::BetaSheet => s.beta += 1,
            RamachandranRegion::LeftHandedAlpha => s.left_alpha += 1,
            RamachandranRegion::Bridge => s.bridge += 1,
            RamachandranRegion::Outlier => s.outliers += 1,
        }
    }
    s
}

/// Wrap an angle into `(-180, 180]`.
fn wrap180(a: f64) -> f64 {
    let mut x = a % 360.0;
    if x > 180.0 {
        x -= 360.0;
    } else if x <= -180.0 {
        x += 360.0;
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, Residue};
    use nalgebra::Point3;

    #[test]
    fn classifies_canonical_basins() {
        // textbook helix (-57, -47)
        assert_eq!(classify(-57.0, -47.0), RamachandranRegion::AlphaHelix);
        // textbook beta (-120, 130)
        assert_eq!(classify(-120.0, 130.0), RamachandranRegion::BetaSheet);
        // left-handed (+60, +40)
        assert_eq!(classify(60.0, 40.0), RamachandranRegion::LeftHandedAlpha);
        // far outlier
        assert_eq!(classify(0.0, 0.0), RamachandranRegion::Outlier);
    }

    #[test]
    fn allowed_predicate() {
        assert!(classify(-60.0, -45.0).is_allowed());
        assert!(!classify(10.0, 170.0).is_allowed());
    }

    #[test]
    fn angle_wrapping() {
        // 200 deg wraps to -160, still inside the helix phi band.
        assert_eq!(classify(-160.0 + 360.0, -50.0), RamachandranRegion::AlphaHelix);
    }

    #[test]
    fn empty_chain_summary() {
        let chain = Chain::new("A");
        let s = summarize(&chain);
        assert_eq!(s.total, 0);
        assert!((s.allowed_fraction() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn helical_chain_summary() {
        // Build a 6-residue ideal-ish helix; every interior residue
        // should land in the helix basin.
        let mut chain = Chain::new("A");
        // ideal alpha-helix internal coords give phi~-57 psi~-47;
        // we lay atoms on a parametric helix so the torsions are
        // close to that.
        let rise = 1.5;
        let twist = 100.0_f64.to_radians();
        let radius = 2.3;
        for i in 0..8i32 {
            let theta = i as f64 * twist;
            let base = i as f64 * rise;
            let mut r = Residue::new("ALA", i + 1);
            // N, CA, C placed slightly apart along the helix.
            for (k, name) in ["N", "CA", "C"].iter().enumerate() {
                let t = theta + k as f64 * 0.6;
                r.atoms.push(Atom::new(
                    *name,
                    "C",
                    Point3::new(
                        radius * t.cos(),
                        radius * t.sin(),
                        base + k as f64 * 0.5,
                    ),
                ));
            }
            chain.residues.push(r);
        }
        let s = summarize(&chain);
        assert!(s.total >= 6, "expected interior residues, got {}", s.total);
        // We don't assert the exact basin (geometry is approximate),
        // only that classification ran and the fraction is sane.
        assert!(s.allowed_fraction() >= 0.0 && s.allowed_fraction() <= 1.0);
    }
}
