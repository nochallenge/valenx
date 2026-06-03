//! A compact backbone-independent rotamer library.
//!
//! A *rotamer* is a discrete, statistically-common sidechain
//! conformation — a tuple of χ (chi) dihedral angles. The classical
//! sidechain-placement algorithms (SCWRL-class repacking, Rosetta
//! `fixbb`-class design) all work by choosing, per residue, one
//! rotamer from a precomputed library and scoring the combination.
//!
//! This module ships a **compact backbone-independent library**: for
//! each amino acid, the rotatable χ angles are enumerated over the
//! three canonical staggered values (gauche+ ≈ +60°, trans ≈ 180°,
//! gauche− ≈ −60°), the choice the Penultimate / Lovell rotamer
//! libraries are built around. A real production pipeline uses the
//! Dunbrack *backbone-dependent* library (rotamer frequencies
//! conditioned on φ/ψ); this compact set is the honest v1 — the same
//! staggered-rotamer idea, without the backbone-conditioned
//! frequencies.
//!
//! From a rotamer the module reconstructs the **sidechain heavy-atom
//! centroid** relative to the backbone — enough for the
//! centroid-resolution packing and design scores. Full all-atom
//! sidechain reconstruction (every Cγ/Cδ/…) is a documented v1
//! simplification: the centroid carries the steric bulk the
//! knowledge-based score needs.

use nalgebra::{Point3, Vector3};
use serde::{Deserialize, Serialize};

use crate::aa::chi_count;
use crate::model::{ideal, ModelResidue};

/// One rotamer: an amino-acid identity, its χ angles (degrees) and a
/// prior probability (the rotamer's relative frequency, used to bias
/// the search toward common conformations).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Rotamer {
    /// One-letter amino-acid code.
    pub aa: char,
    /// χ dihedral angles in degrees, `chi.len() == chi_count(aa)`.
    pub chi: Vec<f64>,
    /// Relative prior probability of this rotamer (`(0, 1]`).
    pub probability: f64,
}

impl Rotamer {
    /// The canonical staggered χ values nearest a rotamer index:
    /// `0 → +60°`, `1 → 180°`, `2 → −60°`.
    fn staggered(index: usize) -> f64 {
        match index % 3 {
            0 => 60.0,
            1 => 180.0,
            _ => -60.0,
        }
    }
}

/// The full rotamer set for one amino-acid identity.
///
/// For an amino acid with `k` χ angles this is the `3^k` combinations
/// of staggered values (capped — see [`rotamers_for`]). Alanine and
/// glycine have a single trivial "rotamer" (no χ).
pub fn rotamers_for(aa: char) -> Vec<Rotamer> {
    let k = chi_count(aa);
    if k == 0 {
        return vec![Rotamer {
            aa,
            chi: Vec::new(),
            probability: 1.0,
        }];
    }
    // Cap the combinatorial blow-up: only the first 2 χ angles get
    // the full 3-way staggered choice, deeper χ angles are fixed at
    // trans. A real library would carry all combinations with
    // measured frequencies; this keeps the per-residue set small
    // (≤ 9) while covering the dominant χ1/χ2 variation.
    let free = k.min(2);
    let combos = 3usize.pow(free as u32);
    let mut out = Vec::with_capacity(combos);
    for combo in 0..combos {
        let mut chi = Vec::with_capacity(k);
        let mut c = combo;
        for _ in 0..free {
            chi.push(Rotamer::staggered(c % 3));
            c /= 3;
        }
        // Deeper χ angles are fixed at trans (180°).
        chi.extend(std::iter::repeat_n(180.0, k - free));
        // Trans χ1 is the most common; weight it up slightly.
        let probability = if chi.first().map(|&x| x == 180.0).unwrap_or(true) {
            1.0
        } else {
            0.7
        };
        out.push(Rotamer {
            aa,
            chi,
            probability,
        });
    }
    // Normalise the priors so they sum to 1.
    let total: f64 = out.iter().map(|r| r.probability).sum();
    if total > 0.0 {
        for r in &mut out {
            r.probability /= total;
        }
    }
    out
}

/// Reconstructs the β-carbon position from the backbone `N`, `CA`,
/// `C` atoms with idealised tetrahedral geometry. This is the
/// standard rebuilt-Cβ formula (Cβ on the bisector of the N–CA–C
/// angle, tilted out of the backbone plane), correct to within the
/// small spread of real Cβ positions.
pub fn rebuild_cb(n: &Point3<f64>, ca: &Point3<f64>, c: &Point3<f64>) -> Point3<f64> {
    let ca_n = (n - ca).normalize();
    let ca_c = (c - ca).normalize();
    // In-plane bisector pointing away from N and C.
    let bisector = -(ca_n + ca_c).normalize();
    // Out-of-plane normal.
    let normal = ca_n.cross(&ca_c).normalize();
    // Cβ sits ~54.75° between the bisector and the normal — the
    // tetrahedral tilt.
    let dir = (bisector * 0.5774 + normal * 0.8165).normalize();
    ca + dir * ideal::CA_CB
}

/// The sidechain-centroid offset, in ångström, for a residue's
/// amino-acid identity and chosen rotamer.
///
/// Returns the distance the sidechain heavy-atom centroid sits
/// *beyond* the Cβ along the Cα→Cβ direction. Larger sidechains push
/// the centroid further out; the χ angles tilt it, modelled here as a
/// small lateral shift. This is the centroid-resolution sidechain the
/// knowledge-based score consumes — see the module note on the
/// all-atom simplification.
fn centroid_reach(aa: char) -> f64 {
    // Roughly: cube-root of the sidechain volume sets the reach.
    let v = crate::aa::sidechain_volume(aa);
    if v <= 0.0 {
        0.0
    } else {
        0.6 + 0.18 * v.cbrt()
    }
}

/// Places the sidechain-centroid pseudo-atom for a residue given its
/// chosen rotamer.
///
/// Builds (or reuses) the Cβ, then offsets the centroid along the
/// Cα→Cβ axis by an amino-acid-specific reach, with a lateral tilt set by the
/// rotamer's first χ angle. Glycine returns its Cα (no sidechain).
///
/// Returns `None` if the residue lacks the `N`/`CA`/`C` backbone
/// needed to anchor the sidechain.
pub fn place_sidechain_centroid(residue: &ModelResidue, rotamer: &Rotamer) -> Option<Point3<f64>> {
    let (n, ca, c) = (residue.n?, residue.ca?, residue.c?);
    if residue.aa == 'G' {
        return Some(ca);
    }
    let cb = residue.cb.unwrap_or_else(|| rebuild_cb(&n, &ca, &c));
    let axis = (cb - ca).normalize();
    let reach = centroid_reach(residue.aa);
    // A lateral tilt direction perpendicular to the Cα→Cβ axis.
    let mut perp = axis.cross(&Vector3::new(0.0, 0.0, 1.0));
    if perp.norm() < 1e-6 {
        perp = axis.cross(&Vector3::new(0.0, 1.0, 0.0));
    }
    perp.normalize_mut();
    let chi1 = rotamer.chi.first().copied().unwrap_or(180.0).to_radians();
    let tilt = 0.35 * reach;
    Some(cb + axis * reach + perp * (tilt * chi1.sin()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alanine_has_one_trivial_rotamer() {
        let r = rotamers_for('A');
        assert_eq!(r.len(), 1);
        assert!(r[0].chi.is_empty());
    }

    #[test]
    fn leucine_rotamer_set_is_bounded() {
        // Leu has 2 χ → 9 staggered combinations.
        let r = rotamers_for('L');
        assert_eq!(r.len(), 9);
        for rot in &r {
            assert_eq!(rot.chi.len(), 2);
        }
        let total: f64 = r.iter().map(|x| x.probability).sum();
        assert!((total - 1.0).abs() < 1e-9, "priors normalised");
    }

    #[test]
    fn lysine_caps_deep_chi_at_trans() {
        // Lys has 4 χ; only first 2 vary → 9 rotamers, χ3/χ4 = 180.
        let r = rotamers_for('K');
        assert_eq!(r.len(), 9);
        for rot in &r {
            assert_eq!(rot.chi.len(), 4);
            assert_eq!(rot.chi[2], 180.0);
            assert_eq!(rot.chi[3], 180.0);
        }
    }

    #[test]
    fn rebuilt_cb_has_right_bond_length() {
        let n = Point3::new(0.0, 1.45, 0.0);
        let ca = Point3::new(0.0, 0.0, 0.0);
        let c = Point3::new(1.52, -0.5, 0.0);
        let cb = rebuild_cb(&n, &ca, &c);
        let len = (cb - ca).norm();
        assert!((len - ideal::CA_CB).abs() < 1e-9, "Cβ bond {len}");
    }

    #[test]
    fn glycine_centroid_is_its_ca() {
        let mut res = ModelResidue::empty('G');
        res.n = Some(Point3::new(0.0, 1.45, 0.0));
        res.ca = Some(Point3::new(0.0, 0.0, 0.0));
        res.c = Some(Point3::new(1.52, -0.5, 0.0));
        let rot = &rotamers_for('G')[0];
        let cen = place_sidechain_centroid(&res, rot).expect("centroid");
        assert_eq!(cen, res.ca.unwrap());
    }
}
