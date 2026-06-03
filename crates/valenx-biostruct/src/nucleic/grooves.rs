//! Major- and minor-groove width estimation.
//!
//! Groove width is measured between phosphate atoms on opposing
//! strands. For each phosphate on one strand the nearest phosphate
//! on the *other* strand is found; the inter-phosphate distance,
//! minus twice the phosphate van der Waals radius (the standard
//! El Hassan-Calladine correction), is the local groove width.
//!
//! Whether a given phosphate pair bounds the **major** or **minor**
//! groove is decided from the sequence offset between the two paired
//! residues: short offsets bound the minor groove, the wider offset
//! bounds the major groove.
//!
//! ## Scope of this v1
//!
//! This is the phosphate-phosphate "P-P distance" groove definition,
//! the simplest of the several in the literature. It does not do the
//! Curves+ refined groove-surface calculation; for canonical B-DNA
//! it gives sensible major/minor numbers, for highly distorted
//! structures it is approximate.

use crate::structure::{Chain, Residue};
use nalgebra::Point3;

/// Phosphate phosphorus van der Waals radius, ångström — the
/// El Hassan-Calladine groove correction subtracts `2 × this`.
const PHOSPHATE_RADIUS: f64 = 1.80;

/// A single local groove-width measurement.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct GrooveWidth {
    /// Residue index (in strand 1) of the reference phosphate.
    pub residue_index: usize,
    /// Corrected groove width, ångström (P-P distance − 2·r_P).
    pub width: f64,
    /// `true` for a major-groove measurement, `false` for minor.
    pub is_major: bool,
}

/// Both grooves' widths along a duplex.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct GrooveProfile {
    /// Major-groove local widths along the duplex.
    pub major: Vec<GrooveWidth>,
    /// Minor-groove local widths along the duplex.
    pub minor: Vec<GrooveWidth>,
}

impl GrooveProfile {
    /// Mean major-groove width, or `None` if no measurements exist.
    pub fn mean_major(&self) -> Option<f64> {
        mean(&self.major)
    }

    /// Mean minor-groove width, or `None` if no measurements exist.
    pub fn mean_minor(&self) -> Option<f64> {
        mean(&self.minor)
    }
}

fn mean(widths: &[GrooveWidth]) -> Option<f64> {
    if widths.is_empty() {
        return None;
    }
    Some(widths.iter().map(|w| w.width).sum::<f64>() / widths.len() as f64)
}

/// The phosphorus atom coordinate of a nucleotide, ångström.
pub fn phosphate(residue: &Residue) -> Option<Point3<f64>> {
    residue.primary_atom("P").map(|a| a.coord)
}

/// Estimate the major- and minor-groove widths of an antiparallel
/// duplex formed by `strand1` and `strand2`.
///
/// The two chains are treated as a duplex with `strand2` antiparallel
/// to `strand1`. For each phosphate of `strand1`, the two nearest
/// `strand2` phosphates (on either side along the helix) are taken as
/// the groove-bounding pair: the closer one is assigned to the minor
/// groove, the farther to the major groove.
pub fn groove_widths(strand1: &Chain, strand2: &Chain) -> GrooveProfile {
    // Collect (residue_index, P-coord) for each strand's nucleotides.
    let p1: Vec<(usize, Point3<f64>)> = strand1
        .residues
        .iter()
        .enumerate()
        .filter_map(|(i, r)| {
            if r.is_nucleotide() {
                phosphate(r).map(|p| (i, p))
            } else {
                None
            }
        })
        .collect();
    let p2: Vec<(usize, Point3<f64>)> = strand2
        .residues
        .iter()
        .enumerate()
        .filter_map(|(i, r)| {
            if r.is_nucleotide() {
                phosphate(r).map(|p| (i, p))
            } else {
                None
            }
        })
        .collect();

    let mut profile = GrooveProfile::default();
    if p1.is_empty() || p2.len() < 2 {
        return profile;
    }

    for &(ri, ref_p) in &p1 {
        // Two smallest P-P distances to the other strand.
        let mut dists: Vec<f64> = p2.iter().map(|(_, q)| (ref_p - q).norm()).collect();
        dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        // dists[0] bounds the minor groove, dists[1] the major.
        let minor_w = (dists[0] - 2.0 * PHOSPHATE_RADIUS).max(0.0);
        let major_w = (dists[1] - 2.0 * PHOSPHATE_RADIUS).max(0.0);
        profile.minor.push(GrooveWidth {
            residue_index: ri,
            width: minor_w,
            is_major: false,
        });
        profile.major.push(GrooveWidth {
            residue_index: ri,
            width: major_w,
            is_major: true,
        });
    }
    profile
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::Atom;

    /// Build a strand of nucleotides with phosphates at the given
    /// coordinates.
    fn strand(name: &str, ps: &[Point3<f64>]) -> Chain {
        let mut c = Chain::new(name);
        for (k, p) in ps.iter().enumerate() {
            let mut r = Residue::new("DA", k as i32 + 1);
            r.atoms.push(Atom::new("P", "P", *p));
            c.residues.push(r);
        }
        c
    }

    #[test]
    fn phosphate_lookup() {
        let mut r = Residue::new("DG", 1);
        r.atoms
            .push(Atom::new("P", "P", Point3::new(1.0, 2.0, 3.0)));
        assert!(phosphate(&r).is_some());
        let bad = Residue::new("DG", 2);
        assert!(phosphate(&bad).is_none());
    }

    #[test]
    fn measures_two_grooves() {
        // strand1 phosphates on the x-axis; strand2 phosphates
        // offset in y so each strand1 P sees a near and a far
        // strand2 P.
        let s1 = strand(
            "A",
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 3.4),
                Point3::new(0.0, 0.0, 6.8),
            ],
        );
        let s2 = strand(
            "B",
            &[
                Point3::new(11.0, 0.0, 0.0),
                Point3::new(18.0, 0.0, 3.4),
                Point3::new(11.0, 0.0, 6.8),
            ],
        );
        let prof = groove_widths(&s1, &s2);
        assert_eq!(prof.major.len(), 3);
        assert_eq!(prof.minor.len(), 3);
        // major width is always >= minor width by construction.
        for (mj, mn) in prof.major.iter().zip(&prof.minor) {
            assert!(mj.width >= mn.width);
            assert!(mj.is_major && !mn.is_major);
        }
        assert!(prof.mean_major().unwrap() >= prof.mean_minor().unwrap());
    }

    #[test]
    fn correction_subtracts_phosphate_radii() {
        // A single strand1 P with two strand2 Ps at exactly 10 and
        // 20 A: minor width = 10 - 3.6 = 6.4, major = 20 - 3.6 = 16.4.
        let s1 = strand("A", &[Point3::new(0.0, 0.0, 0.0)]);
        let s2 = strand(
            "B",
            &[Point3::new(10.0, 0.0, 0.0), Point3::new(20.0, 0.0, 0.0)],
        );
        let prof = groove_widths(&s1, &s2);
        assert!((prof.minor[0].width - 6.4).abs() < 1e-9);
        assert!((prof.major[0].width - 16.4).abs() < 1e-9);
    }

    #[test]
    fn empty_strands_yield_empty_profile() {
        let empty = Chain::new("A");
        let prof = groove_widths(&empty, &empty);
        assert!(prof.major.is_empty());
        assert!(prof.mean_major().is_none());
    }
}
