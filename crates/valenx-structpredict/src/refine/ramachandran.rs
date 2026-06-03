//! **Feature 15 — backbone-dihedral / Ramachandran refinement.**
//!
//! A model's backbone is only physically plausible if its φ/ψ
//! dihedral angles fall in the **allowed regions of the Ramachandran
//! plot** — the basins (right-handed α, β, left-handed α,
//! polyproline) that real protein backbones actually occupy. A
//! freshly-built or fragment-assembled model often has a few residues
//! in disallowed regions ("Ramachandran outliers").
//!
//! This module refines them: it computes each residue's φ/ψ, finds
//! the residues sitting in disallowed space, and **snaps each outlier
//! to the nearest allowed basin**, then rebuilds the affected
//! backbone from the corrected torsions. The result has every residue
//! in (or near) an allowed region.
//!
//! Allowed-region membership uses [`valenx_biostruct`]'s Ramachandran
//! classifier where a [`valenx_biostruct::Structure`] is available;
//! this module additionally provides a coordinate-free basin snap so
//! it can refine a [`crate::model::ProteinModel`] directly.

use serde::{Deserialize, Serialize};

use crate::error::{Result, StructPredictError};
use crate::model::{build_backbone_from_torsions, ProteinModel};

/// The outcome of a Ramachandran refinement pass.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RamachandranRefinement {
    /// Residue indices that were Ramachandran outliers before
    /// refinement.
    pub outliers_before: Vec<usize>,
    /// Residue indices still outliers after refinement (basins very
    /// near a disallowed edge can remain).
    pub outliers_after: Vec<usize>,
    /// φ/ψ angles (degrees) that were snapped, as
    /// `(residue, old_phi, old_psi, new_phi, new_psi)`.
    pub adjustments: Vec<(usize, f64, f64, f64, f64)>,
}

impl RamachandranRefinement {
    /// Number of outliers removed.
    pub fn removed(&self) -> usize {
        self.outliers_before.len().saturating_sub(self.outliers_after.len())
    }
}

/// The canonical allowed (φ, ψ) basin centres, degrees.
const ALLOWED_BASINS: &[(f64, f64)] = &[
    (-63.0, -42.0),  // right-handed α-helix
    (-120.0, 130.0), // β-strand
    (-75.0, 145.0),  // polyproline-II
    (57.0, 47.0),    // left-handed α
    (-90.0, 0.0),    // bridge / turn region
];

/// Whether a (φ, ψ) pair (degrees) is in an allowed Ramachandran
/// region — within a generous radius of one of the canonical basins.
///
/// The radius (`~60°`) is deliberately generous: it accepts the
/// broad allowed regions of a real Ramachandran plot and flags only
/// genuinely disallowed geometry.
pub fn is_allowed(phi: f64, psi: f64) -> bool {
    ALLOWED_BASINS
        .iter()
        .any(|&(bp, bq)| angular_distance(phi, bp).hypot(angular_distance(psi, bq)) < 60.0)
}

/// The smallest signed angular difference `a − b`, wrapped to
/// `(−180°, 180°]`.
fn angular_distance(a: f64, b: f64) -> f64 {
    let mut d = (a - b) % 360.0;
    if d > 180.0 {
        d -= 360.0;
    } else if d <= -180.0 {
        d += 360.0;
    }
    d
}

/// The allowed basin nearest a (φ, ψ) pair.
fn nearest_basin(phi: f64, psi: f64) -> (f64, f64) {
    let mut best = ALLOWED_BASINS[0];
    let mut best_d = f64::INFINITY;
    for &(bp, bq) in ALLOWED_BASINS {
        let d = angular_distance(phi, bp).hypot(angular_distance(psi, bq));
        if d < best_d {
            best_d = d;
            best = (bp, bq);
        }
    }
    best
}

/// Computes the per-residue (φ, ψ) backbone dihedrals of a model.
///
/// Residue `0` has no φ (no preceding C) and the last residue no ψ
/// (no following N); those are reported as `0.0`. A residue without a
/// full backbone yields `(0.0, 0.0)`.
pub fn model_phi_psi(model: &ProteinModel) -> Vec<(f64, f64)> {
    let n = model.residues.len();
    let mut out = vec![(0.0, 0.0); n];
    for (i, slot) in out.iter_mut().enumerate() {
        let res = &model.residues[i];
        let (Some(ni), Some(cai), Some(ci)) = (res.n, res.ca, res.c) else {
            continue;
        };
        // φ = dihedral C(i-1) - N(i) - CA(i) - C(i).
        let phi = if i > 0 {
            if let Some(prev_c) = model.residues[i - 1].c {
                dihedral_deg(prev_c, ni, cai, ci)
            } else {
                0.0
            }
        } else {
            0.0
        };
        // ψ = dihedral N(i) - CA(i) - C(i) - N(i+1).
        let psi = if i + 1 < n {
            if let Some(next_n) = model.residues[i + 1].n {
                dihedral_deg(ni, cai, ci, next_n)
            } else {
                0.0
            }
        } else {
            0.0
        };
        *slot = (phi, psi);
    }
    out
}

/// The signed dihedral angle (degrees) of four points.
fn dihedral_deg(
    a: nalgebra::Point3<f64>,
    b: nalgebra::Point3<f64>,
    c: nalgebra::Point3<f64>,
    d: nalgebra::Point3<f64>,
) -> f64 {
    let b1 = b - a;
    let b2 = c - b;
    let b3 = d - c;
    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    let m = n1.cross(&b2.normalize());
    let x = n1.dot(&n2);
    let y = m.dot(&n2);
    y.atan2(x).to_degrees()
}

/// Refines a model's backbone by snapping Ramachandran outliers to
/// the nearest allowed basin.
///
/// Computes every residue's φ/ψ, identifies the outliers, snaps each
/// outlier's torsions to the nearest allowed basin, and rebuilds the
/// whole backbone from the corrected torsion set. The first and last
/// residues' missing φ/ψ are left at a default helical value when
/// rebuilding.
///
/// # Errors
/// [`StructPredictError::Invalid`] for a model without a complete
/// backbone (φ/ψ cannot be measured otherwise).
pub fn refine_ramachandran(model: &mut ProteinModel) -> Result<RamachandranRefinement> {
    if !model.is_complete() {
        return Err(StructPredictError::invalid(
            "model",
            "Ramachandran refinement needs a complete backbone",
        ));
    }
    let phi_psi = model_phi_psi(model);
    let n = phi_psi.len();
    if n < 3 {
        return Err(StructPredictError::invalid(
            "model",
            "need at least 3 residues",
        ));
    }

    // Identify outliers (skip the first/last, whose φ/ψ are partial).
    let mut outliers_before = Vec::new();
    for (i, &(phi, psi)) in phi_psi.iter().enumerate() {
        if i == 0 || i == n - 1 {
            continue;
        }
        if !is_allowed(phi, psi) {
            outliers_before.push(i);
        }
    }

    // Build a corrected torsion set: outliers snap to the nearest
    // basin; everyone else keeps their measured angles.
    let mut torsions: Vec<(f64, f64)> = phi_psi.clone();
    let mut adjustments = Vec::new();
    for &i in &outliers_before {
        let (old_phi, old_psi) = phi_psi[i];
        let (new_phi, new_psi) = nearest_basin(old_phi, old_psi);
        torsions[i] = (new_phi, new_psi);
        adjustments.push((i, old_phi, old_psi, new_phi, new_psi));
    }
    // Give the terminal residues a sensible helical default so the
    // rebuild is well-defined.
    if n > 0 {
        torsions[0] = (-63.0, torsions[0].1);
        torsions[n - 1] = (torsions[n - 1].0, -42.0);
    }

    build_backbone_from_torsions(model, &torsions)?;

    // Recheck.
    let after_pp = model_phi_psi(model);
    let mut outliers_after = Vec::new();
    for (i, &(phi, psi)) in after_pp.iter().enumerate() {
        if i == 0 || i == n - 1 {
            continue;
        }
        if !is_allowed(phi, psi) {
            outliers_after.push(i);
        }
    }

    Ok(RamachandranRefinement {
        outliers_before,
        outliers_after,
        adjustments,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helix_is_all_allowed() {
        // α-helix (φ, ψ) ≈ (-63, -42) is firmly allowed.
        assert!(is_allowed(-63.0, -42.0));
        assert!(is_allowed(-120.0, 130.0)); // β-strand
                                            // A clearly disallowed point.
        assert!(!is_allowed(0.0, 90.0));
    }

    #[test]
    fn refinement_removes_outliers() {
        // Build a chain that is mostly helical but has a disallowed
        // residue spliced in.
        let mut m = ProteinModel::from_sequence("AAAAAAAAAA").expect("model");
        let mut torsions = vec![(-63.0, -42.0); 10];
        torsions[5] = (10.0, 90.0); // disallowed
        build_backbone_from_torsions(&mut m, &torsions).expect("build");
        let refinement = refine_ramachandran(&mut m).expect("refine");
        assert!(
            !refinement.outliers_before.is_empty(),
            "an outlier was introduced"
        );
        assert!(
            refinement.outliers_after.len() < refinement.outliers_before.len(),
            "refinement removed outliers"
        );
        assert!(refinement.removed() >= 1);
    }

    #[test]
    fn angular_distance_wraps() {
        // 170° and -170° are only 20° apart, not 340°.
        assert!((angular_distance(170.0, -170.0).abs() - 20.0).abs() < 1e-9);
    }

    #[test]
    fn incomplete_model_rejected() {
        let m = ProteinModel::from_sequence("AAAA").expect("model"); // no coords
        let mut m = m;
        assert!(refine_ramachandran(&mut m).is_err());
    }
}
